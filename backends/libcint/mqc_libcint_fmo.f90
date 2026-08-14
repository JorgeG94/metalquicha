!! The fragment molecular orbital method, two-body, for non-covalent fragments
module mqc_libcint_fmo
   !! FMO2: two nested self-consistencies, then a pass over pairs.
   !!
   !! **The inner SCF** is an ordinary fragment SCF -- orbitals converged
   !! against a fixed external potential.
   !!
   !! **The outer SCF**, the monomer loop, converges that potential. Every
   !! fragment is solved in the field of all the others, which changes its
   !! density, which changes the field every other fragment sees. Iterating the
   !! two until the monomer energies stop moving makes the fragment densities
   !! mutually consistent, and it is the step that separates this from computing
   !! fragments independently and adding corrections afterwards.
   !!
   !! **What the field is.** The embedding operator on fragment X, added to its
   !! core Hamiltonian and so present in every Fock build of its inner SCF:
   !!
   !!     u^X_mn = sum_{K/=X} [ -sum_{A in K} Z_A <m|1/|r-R_A||n>
   !!                           + sum_{ls in K} D^K_ls (mn|ls) ]
   !!
   !! -- exact nuclear attraction to the other fragments' nuclei, and the exact
   !! Coulomb operator of their electron densities. Both are integrals over the
   !! real densities, which is the point: a fragment feels the field the rest of
   !! the system genuinely makes, and a dimer feels the field of the converged
   !! monomers, not a set of fitted point charges.
   !!
   !! Replacing that second term with atomic point charges is the ESP-PTC
   !! approximation, available as `esp = "ptc"`. In FMO that is reserved for
   !! *distant* fragments.
   !!
   !! **Three methods live here, along two independent axes.** `esp` decides
   !! what field the fragments were computed in; `expansion` decides how their
   !! energies are assembled. The named methods are particular pairings:
   !!
   !!     esp = "exact", expansion = "fmo"   FMO2
   !!     esp = "ptc",   expansion = "mbe"   electrostatically embedded MBE
   !!     esp = "none",  expansion = "mbe"   plain MBE
   !!
   !! Keeping the axes apart matters because point charges alone do not make a
   !! calculation EE-MBE. ESP-PTC under the FMO expansion is still FMO -- it
   !! carries internal energies and the response term -- and it is a genuinely
   !! different quantity, by
   !! `Tr(D_I u_I) + Tr(D_J u_J) - Tr((D_I+D_J) u_IJ)`, which does not vanish.
   !! On water trimers the two land within about 15% of each other's error, so
   !! nothing loud goes wrong if they are confused; they are simply not the same
   !! method, and `expansion` exists so that either can be asked for by name.
   !! See `validation/sweep_fmo`.
   !!
   !! **The energy.** With `E'` an internal energy -- the fragment's own energy
   !! with its polarised density, not counting its interaction with the field:
   !!
   !!     E = sum_I E'_I + sum_{I<J} (E'_IJ - E'_I - E'_J)
   !!                    + sum_{I<J} Tr(dD_IJ u_IJ)
   !!
   !! The pair term carries the whole I-J interaction, since the dimer is solved
   !! with both present. The last is the pair's density response to the field of
   !! everything outside it, `dD_IJ` being the dimer density less the two
   !! monomer densities laid side by side.
   !!
   !! **Non-covalent fragments only, and this is enforced.** Every fragment must
   !! be a whole molecule. None of the machinery for cutting a bond -- adjusted
   !! fragment orbitals, hybrid orbital projection, hydrogen caps -- is present,
   !! so a partition that severs one is refused rather than answered.
   !!
   !! The check is not a formality. Cutting a single bond leaves both fragments
   !! with an odd electron count, which the closed-shell check would catch on
   !! its own; but cut an even number per fragment -- a ring, a double bond --
   !! and every count stays even. Cyclopropane split into three CH2 used to come
   !! back 0.28 Hartree low, which is 176 kcal/mol in the shape of an answer.
   !! Connectivity is therefore checked directly, with the criterion
   !! [[mqc_bond_perception]] uses everywhere else.
   !!
   !! **Why two fragments is the test that matters.** With two fragments nothing
   !! lies outside the dimer, so `u_12` vanishes, the monomer terms cancel
   !! algebraically, and `E = E'_12` -- an ordinary supermolecular SCF. FMO2 on
   !! two fragments is therefore not an approximation at all, and disagreement
   !! with a plain RHF is a bug in how fragments are built rather than an error
   !! in the method. See `validation/check_fmo`.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use pic_mpi_lib, only: comm_t, allreduce, MPI_SUM
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_elements, only: element_vdw_radius
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_bond_perception, only: connected_components
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_direct, only: schwarz_bounds, build_fock_direct, direct_stats_t
   use mqc_libcint_esp, only: esp_matrices
   use mqc_libcint_charges, only: mulliken_charges, chelpg_charges
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   implicit none
   private

   public :: fmo_options_t
   public :: fmo_result_t
   public :: run_fmo2

   type :: fmo_options_t
      !! What to run, and how hard
      character(len=64) :: basis = "6-31g"
      character(len=16) :: esp = "exact"
         !! How a fragment's neighbours are represented to it.
         !!
         !! `"exact"` builds the embedding from the neighbours' actual
         !! densities: nuclear attraction integrals plus the Coulomb operator of
         !! their density matrices. This is FMO.
         !!
         !! `"ptc"` collapses each neighbour to atomic point charges
         !! (`ptc_charges` picks which). In FMO this approximation is reserved
         !! for distant fragments; using it everywhere is electrostatically
         !! embedded MBE, which is a different and perfectly good method.
         !!
         !! `"none"` removes the embedding, leaving a plain many-body expansion
         !! -- the baseline that says what embedding buys.
         !!
         !! **Measured** against supermolecular RHF on stacked water clusters
         !! (`validation/sweep_fmo`), 5 monomers, error in Hartree:
         !!
         !!     O-O sep    exact ESP    ptc-mulliken   ptc-chelpg    no embed
         !!     2.70 A     -4.9e-04       7.1e-05      -2.1e-04      3.7e-03
         !!     2.90 A     -1.1e-04       2.0e-04      -3.1e-05      2.9e-03
         !!     3.20 A     -1.0e-06       1.6e-04       1.3e-05      1.9e-03
         !!     4.00 A      8.6e-07       2.6e-05       1.1e-06      4.6e-04
         !!     6.00 A      1.0e-12       9.1e-07       8.1e-10      3.7e-05
         !!     9.00 A      1.7e-13       8.1e-08       2.0e-10      3.4e-06
         !!
         !! Read the last two rows first. Once the fragments are far enough
         !! apart that exchange and charge transfer between them have died --
         !! while the electrostatic field very much has not -- FMO2 with the
         !! exact ESP is *exact*, to 1e-13. Nothing but a correct embedding
         !! operator does that. The point-charge variants plateau five orders of
         !! magnitude short and stay there however far the fragments are moved,
         !! because what is left is the charge approximation itself.
         !!
         !! The short-separation rows read the other way, and the temptation is
         !! to conclude point charges are better there. They are not. The exact
         !! ESP's error at 2.70 A is the genuine three-body term, which FMO2
         !! does not contain and no embedding can supply; it is negative at
         !! every separation and shrinks monotonically as that term dies. The
         !! point-charge error has the opposite sign over part of the range, so
         !! it partly cancels against the three-body term, and near contact the
         !! cancellation flatters it. That is luck, not accuracy -- and it is
         !! not something to rely on, because the two errors have no reason to
         !! stay matched for a different system.
      character(len=16) :: far_field = "mulliken"
         !! What a distant fragment contributes to the embedding: the level of
         !! theory the long-range part is evaluated at.
         !!
         !! `"mulliken"` -- point charges from Mulliken populations. The
         !! default, and what production FMO codes use. The approximation being
         !! made is specifically a population one: the distant fragment's
         !! `sum_ls D_ls (mn|ls)` is replaced by its atomic populations, and
         !! Mulliken populations are what that term reduces to.
         !!
         !! `"chelpg"` -- point charges from an ESP fit. Production FMO codes
         !! do not use these and the choice is not a studied one, so what
         !! follows is a measurement rather than a citation.
         !!
         !! At long separation, where the point-charge approximation is the
         !! *only* error left -- no three-body term to cancel against -- CHELPG
         !! is better than Mulliken by one to two orders of magnitude
         !! (`validation/sweep_fmo`, error against supermolecular RHF, Hartree):
         !!
         !!     3 waters at 9 A     mulliken 2.50e-08    chelpg -5.53e-10
         !!     5 waters at 9 A     mulliken 1.05e-07    chelpg  1.33e-08
         !!
         !! That regime is exactly the one the far field operates in, so the
         !! signal is on point and in the direction the charges themselves
         !! predict: CHELPG is fitted to reproduce a potential, and reproducing
         !! a potential is the job. Near contact the comparison is muddied by
         !! cancellation against the three-body term and should not be read.
         !!
         !! Measured on stacked water clusters in one basis. Promising, not
         !! established -- if it matters to a result, check it there.
         !!
         !! `"ignore"` -- distant fragments contribute nothing at all, nuclei
         !! included. Not an approximation to the field so much as a decision
         !! not to have one past the cutoff, which makes it the honest way to
         !! ask what the long-range field is worth: set the cutoff where you
         !! mean to and compare. It is also the cheapest, since no charges need
         !! computing and no integrals evaluating for the distant atoms.
         !!
         !! Read whenever any fragment is distant -- so with `esp = "ptc"`,
         !! which makes all of them distant, and with `esp = "exact"` for those
         !! beyond `resppc`.
      real(dp) :: resppc = 2.0_dp
         !! Where the exact ESP gives way to point charges, as a unitless
         !! separation. Read only when `esp = "exact"`.
         !!
         !! FMO does not evaluate the four-index Coulomb term between every
         !! pair of fragments. Beyond a separation it approximates the distant
         !! fragment's electron density by its atomic populations, which is what
         !! makes the method scale: the exact term is needed only within a
         !! neighbourhood, so the cost per fragment stops growing once the
         !! system is larger than that neighbourhood.
         !!
         !! The separation is measured the way FMO measures it -- the smallest
         !! interatomic distance between the two fragments, divided by the sum
         !! of the two atoms' van der Waals radii:
         !!
         !!     R_IJ = min over A in I, B in J of |R_A - R_B| / (vdw_A + vdw_B)
         !!
         !! so it is a contact distance rather than a centre-of-mass one, and a
         !! large flat fragment does not count as far away merely because its
         !! middle is. The 2.0 default is GAMESS's RESPPC.
         !!
         !! Negative disables the approximation and makes every fragment exact.
         !! Zero makes every fragment distant, which is then identical to
         !! `esp = "ptc"` -- and is asserted to be, in `check_fmo`, since the
         !! two paths reach it by quite different routes.
         !!
         !! **What the default costs**, against the same calculation with the
         !! approximation off (`validation/sweep_fmo`, 5 waters, Hartree):
         !!
         !!     O-O sep    RESPPC = 2      no cutoff
         !!     2.70 A     -4.90e-04      -4.93e-04     nothing is beyond it yet
         !!     2.90 A     -1.08e-04      -1.11e-04
         !!     3.20 A      1.23e-05      -1.03e-06     it starts to bite
         !!     4.00 A      2.51e-06       8.64e-07
         !!     6.00 A     -1.81e-08       1.36e-12     fully engaged
         !!     9.00 A      8.07e-08       1.71e-13
         !!
         !! Read as ratios that looks alarming -- five orders of magnitude at
         !! 9 Angstrom. Read as energies it is the approximation working as
         !! designed: where the cutoff costs most in relative terms, the
         !! absolute error it leaves is 8e-08 Hartree, which is five hundredths
         !! of a millikcal and nothing at all. Where the error is large enough
         !! to matter, near contact, the cutoff has not engaged and costs
         !! nothing. Precision is being given up precisely where precision is
         !! not the binding constraint, which is the whole bargain.
      character(len=16) :: expansion = "fmo"
         !! How the fragment energies are assembled into a total. Orthogonal to
         !! `esp`, which decides what field they were computed in, because the
         !! two choices are genuinely independent and the named methods are
         !! particular pairings of them:
         !!
         !!     esp = "exact", expansion = "fmo"   FMO2
         !!     esp = "ptc",   expansion = "mbe"   electrostatically embedded MBE
         !!     esp = "none",  expansion = "mbe"   plain MBE
         !!
         !! `"fmo"` sums internal energies and adds the density-response term,
         !! as in the module header. `"mbe"` is the ordinary many-body
         !! expansion over the *total* embedded energies:
         !!
         !!     E = sum_I E_I + sum_{I<J} (E_IJ - E_I - E_J)
         !!
         !! with no internal-energy subtraction and no response term. The two
         !! differ by Tr(D_I u_I) + Tr(D_J u_J) - Tr((D_I+D_J) u_IJ), which does
         !! not vanish: `u_I` is the field on I from everything including J,
         !! while `u_IJ` excludes both. So these are different methods and not
         !! two spellings of one.
         !!
         !! A caveat that belongs to EE-MBE rather than to this implementation:
         !! in the pair difference, the point-charge I-J interaction is removed
         !! twice -- once with `E_I`, which was embedded in J's charges, and
         !! once with `E_J` -- while the dimer supplies the real one. That is
         !! inherent to the expansion, not a defect here, and it is part of why
         !! FMO carries internal energies instead.
      integer :: max_outer = 50
         !! Cap on the outer (monomer) SCF
      real(dp) :: outer_tol = 1.0e-7_dp
         !! Outer convergence, on the sum of monomer energies in Hartree.
         !! The energy rather than the density because it is what the answer is
         !! made of, and rather than the charges because with `esp = "exact"`
         !! there are no charges in the loop to converge.
      integer :: scf_max_iter = 100
      real(dp) :: scf_energy_tol = 1.0e-9_dp
      real(dp) :: scf_density_tol = 1.0e-7_dp
   end type fmo_options_t

   type :: fmo_result_t
      !! The energy, and enough of the parts to see where it came from
      real(dp) :: energy = 0.0_dp                !! The FMO2 total
      real(dp) :: monomer_sum = 0.0_dp           !! sum_I E'_I
      real(dp) :: pair_sum = 0.0_dp              !! sum of the pair corrections
      real(dp) :: response_sum = 0.0_dp          !! sum of Tr(dD u), the last term
      integer :: outer_iterations = 0            !! passes of the monomer SCF
      real(dp) :: outer_change = 0.0_dp          !! last movement of the monomer sum
      logical :: converged = .false.
      real(dp), allocatable :: monomer_energy(:)     !! E'_I
      real(dp), allocatable :: pair_correction(:, :)  !! the (I,J) term, upper triangle
      real(dp), allocatable :: charges(:)            !! Mulliken, for reporting only
   end type fmo_result_t

   !> One fragment, and its place in the whole
   type :: fragment_t
      !! What persists between passes, and nothing that does not
      !!
      !! No molecule and no Schwarz bounds. Those are heavy, they are rebuilt
      !! from the geometry wherever they are needed, and holding them would mean
      !! materialising every fragment up front -- a serial bottleneck, and once
      !! this is distributed, the very objects that must not cross a wire. What
      !! stays is small: a density, an energy, a few indices.
      integer :: nao = 0                         !! size of its basis, for block layout
      real(dp), allocatable :: density(:, :)
      real(dp), allocatable :: charges(:)        !! computed while its molecule was in hand
      integer, allocatable :: atoms(:)           !! its atoms, as system indices
      integer, allocatable :: z(:)
      character(len=2), allocatable :: sym(:)
      real(dp), allocatable :: xyz(:, :)
      real(dp) :: energy = 0.0_dp                !! internal, E'
      real(dp) :: energy_total = 0.0_dp          !! as the SCF reported it, with the field
      integer :: nelec = 0
      integer, allocatable :: near(:)            !! fragments close enough for the exact term
   end type fragment_t

contains

   subroutine run_fmo2(atomic_numbers, symbols, coordinates, owner, opts, res, error, comm)
      !! Run FMO2 over a system already partitioned into whole molecules
      !!
      !! `owner(i)` is the fragment index of atom `i`, numbered 1..n_frag with
      !! no gaps -- which is what `connected_components` in
      !! [[mqc_bond_perception]] produces. Coordinates are Bohr.
      !!
      !! Two phases, and they are the two the distributed path hands out: the
      !! monomers, iterated to self-consistency, then the pairs once.
      integer, intent(in) :: atomic_numbers(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: coordinates(:, :)
      integer, intent(in) :: owner(:)
      type(fmo_options_t), intent(in) :: opts
      type(fmo_result_t), intent(out) :: res
      type(error_t), intent(inout) :: error
      type(comm_t), intent(in), optional :: comm
         !! Spread the fragment work over these ranks. Every rank runs this same
         !! routine on the same geometry and assembles only what it is given, so
         !! nothing heavy is ever sent -- only the densities and energies that
         !! the next pass genuinely needs, which are small.

      type(fragment_t), allocatable :: frag(:)
      integer :: n_atoms, n_frag, i
      logical :: all_converged

      n_atoms = size(atomic_numbers)
      if (size(owner) /= n_atoms .or. size(coordinates, 2) /= n_atoms) then
         call error%set(ERROR_VALIDATION, "fmo: owner and coordinates must cover "// &
                        "every atom")
         return
      end if

      call build_fragments(atomic_numbers, symbols, coordinates, owner, opts, &
                           frag, n_frag, error)
      if (error%has_error()) return

      allocate (res%monomer_energy(n_frag), source=0.0_dp)
      allocate (res%pair_correction(n_frag, n_frag), source=0.0_dp)
      all_converged = .true.

      call calculate_monomers(frag, n_frag, atomic_numbers, coordinates, opts, res, &
                              all_converged, error, comm)
      if (error%has_error()) return

      do i = 1, n_frag
         if (opts%expansion == "mbe") then
            res%monomer_energy(i) = frag(i)%energy_total
         else
            res%monomer_energy(i) = frag(i)%energy
         end if
      end do
      res%monomer_sum = sum(res%monomer_energy)

      call calculate_polymers(frag, n_frag, atomic_numbers, coordinates, opts, res, &
                              all_converged, error, comm)
      if (error%has_error()) return

      res%energy = res%monomer_sum + res%pair_sum + res%response_sum

      ! Reported, not used. What the fragments look like once the field has
      ! settled is the first thing to look at when a number seems wrong.
      call report_charges(frag, n_frag, n_atoms, res%charges, error)
      if (error%has_error()) return

      res%converged = res%converged .and. all_converged
      if (.not. all_converged) then
         call error%set(ERROR_VALIDATION, "fmo: at least one fragment SCF did not "// &
                        "converge, so the total is not trustworthy")
         return
      end if
   end subroutine run_fmo2

   subroutine open_fragment(frag_z, frag_sym, frag_xyz, opts, mol, bounds, error)
      !! A fragment's molecule and Schwarz bounds, built from its geometry
      !!
      !! Called wherever one is needed and discarded straight after. That is the
      !! opposite of caching them, and deliberately: the basis set behind
      !! `build_libcint_molecule` is itself cached, so what remains is arithmetic
      !! that parallelises, where a cache of molecules is memory that has to be
      !! filled serially and cannot be shared between ranks anyway.
      integer, intent(in) :: frag_z(:)
      character(len=2), intent(in) :: frag_sym(:)
      real(dp), intent(in) :: frag_xyz(:, :)
      type(fmo_options_t), intent(in) :: opts
      type(libcint_molecule_t), intent(out) :: mol
      real(dp), allocatable, intent(out) :: bounds(:, :)
      type(error_t), intent(inout) :: error

      call build_libcint_molecule(frag_z, frag_sym, frag_xyz, trim(opts%basis), mol, error)
      if (error%has_error()) return
      call schwarz_bounds(mol, bounds, error)
   end subroutine open_fragment

   subroutine build_fragments(z, symbols, coords, owner, opts, frag, n_frag, error)
      !! Each fragment, and which of the others it needs the exact term for
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: owner(:)
      type(fmo_options_t), intent(in) :: opts
      type(fragment_t), allocatable, intent(out) :: frag(:)
      integer, intent(out) :: n_frag
      type(error_t), intent(inout) :: error

      integer, allocatable :: count_per(:)
      integer :: n_atoms, i, f

      n_atoms = size(z)
      n_frag = maxval(owner)
      if (minval(owner) < 1) then
         call error%set(ERROR_VALIDATION, "fmo: fragment labels start at 1")
         return
      end if
      if (n_frag < 2) then
         call error%set(ERROR_VALIDATION, "fmo: the system is one fragment, so there "// &
                        "is nothing to expand -- run an ordinary SCF")
         return
      end if

      allocate (count_per(n_frag), source=0)
      do i = 1, n_atoms
         count_per(owner(i)) = count_per(owner(i)) + 1
      end do
      if (any(count_per == 0)) then
         call error%set(ERROR_VALIDATION, "fmo: fragment labels have a gap, so one "// &
                        "fragment has no atoms")
         return
      end if

      call refuse_severed_bonds(z, coords, owner, n_atoms, error)
      if (error%has_error()) return

      allocate (frag(n_frag))
      do f = 1, n_frag
         frag(f)%atoms = pack([(i, i=1, n_atoms)], owner == f)
         frag(f)%z = z(frag(f)%atoms)
         frag(f)%sym = symbols(frag(f)%atoms)
         frag(f)%xyz = coords(:, frag(f)%atoms)
         frag(f)%nelec = sum(frag(f)%z)
         if (mod(frag(f)%nelec, 2) /= 0) then
            call error%set(ERROR_VALIDATION, "fmo: fragment "//to_char(f)//" has an odd "// &
                           "electron count and this is a closed-shell method")
            return
         end if

         ! Built only to learn how many basis functions it has, which the
         ! block layout needs everywhere, then dropped. Every rank can work this
         ! out for itself without asking anyone.
         block
            type(libcint_molecule_t) :: probe
            real(dp), allocatable :: probe_bounds(:, :)
            call open_fragment(frag(f)%z, frag(f)%sym, frag(f)%xyz, opts, probe, &
                               probe_bounds, error)
            if (error%has_error()) return
            frag(f)%nao = probe%nao
         end block
      end do

      ! Geometry fixes the near sets, so they are found once rather than each
      ! outer pass.
      do f = 1, n_frag
         call near_fragments(frag, n_frag, [f], effective_resppc(opts), frag(f)%near, error)
         if (error%has_error()) return
      end do
      call logger%verbose("  fmo: fragment 1 treats "//to_char(size(frag(1)%near))// &
                          " of "//to_char(n_frag - 1)//" neighbours exactly")
   end subroutine build_fragments

   subroutine refuse_severed_bonds(z, coords, owner, n_atoms, error)
      !! Refuse a partition that cuts a covalent bond
      !!
      !! This method has no way to represent a cut bond -- no adjusted fragment
      !! orbitals, no hybrid orbital projection, no caps -- so a fragment with a
      !! dangling valence is not an approximation it makes, it is a question it
      !! cannot be asked.
      !!
      !! Worth checking explicitly rather than trusting, because the failure is
      !! quiet in exactly the cases that matter. Cutting a single bond leaves
      !! both fragments with an odd electron count, and the closed-shell check
      !! catches that by accident. Cut an even number per fragment -- a ring, a
      !! double bond -- and every count stays even, nothing objects, and the
      !! answer is wrong by a large margin: cyclopropane split into three CH2
      !! comes back 0.28 Hartree low, which is 176 kcal/mol and looks like a
      !! number rather than a mistake.
      !!
      !! The criterion is the one [[mqc_bond_perception]] uses everywhere else,
      !! reused rather than restated so the two cannot drift.
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: owner(:)
      integer, intent(in) :: n_atoms
      type(error_t), intent(inout) :: error

      type(system_geometry_t) :: geom
      integer, allocatable :: component(:)
      integer :: n_components, i, j

      geom%total_atoms = n_atoms
      allocate (geom%element_numbers(n_atoms), source=z)
      allocate (geom%coordinates(3, n_atoms), source=coords)
      call connected_components(geom, component, n_components)

      ! A covalently connected group spanning two fragments is a severed bond.
      ! Reported by naming one offending pair rather than counting them: the
      ! first one is enough to see what was done wrong.
      do i = 1, n_atoms
         do j = i + 1, n_atoms
            if (component(i) /= component(j)) cycle
            if (owner(i) == owner(j)) cycle
            call error%set(ERROR_VALIDATION, "fmo: the partition cuts a covalent "// &
                           "molecule -- atoms "//to_char(i)//" and "//to_char(j)// &
                           " are covalently connected but were put in fragments "// &
                           to_char(owner(i))//" and "//to_char(owner(j))//". This "// &
                           "method has no way to cap a cut bond, so it cannot answer "// &
                           "for that partition; fragment on whole molecules")
            return
         end do
      end do
   end subroutine refuse_severed_bonds

   function effective_resppc(opts) result(r)
      !! The cutoff actually in force
      !!
      !! With `esp = "ptc"` every fragment is distant by construction, which is
      !! a cutoff of zero -- so the two ways of asking for an all-point-charge
      !! embedding meet here rather than in two separate code paths that could
      !! drift apart.
      type(fmo_options_t), intent(in) :: opts
      real(dp) :: r

      r = opts%resppc
      if (opts%esp == "ptc") r = 0.0_dp
      if (opts%esp == "none") r = 0.0_dp
   end function effective_resppc

   subroutine all_charges(frag, n_frag, n_atoms, opts, q_all, error)
      !! One partial charge per atom, from the fragment densities on hand
      !!
      !! Needed whenever any fragment is distant enough to be approximated,
      !! which under the default `resppc` is most of them in a large system.
      type(fragment_t), intent(in) :: frag(:)
      integer, intent(in) :: n_frag, n_atoms
      type(fmo_options_t), intent(in) :: opts
      real(dp), allocatable, intent(out) :: q_all(:)
      type(error_t), intent(inout) :: error

      integer :: f

      allocate (q_all(n_atoms), source=0.0_dp)
      if (opts%esp == "none") return
      ! Nothing reads these if distant fragments contribute nothing, and for
      ! CHELPG they are not cheap enough to compute on the off chance.
      if (opts%far_field == "ignore") return
      ! Taken from what each fragment recorded while its molecule existed, not
      ! recomputed -- rebuilding every fragment just to re-derive charges it
      ! already knows would undo the point of building on demand.
      do f = 1, n_frag
         if (.not. allocated(frag(f)%charges)) cycle
         q_all(frag(f)%atoms) = frag(f)%charges
      end do
   end subroutine all_charges

   subroutine embedding_operator(mol, group_z, group_sym, group_xyz, near, frag, &
                                 n_frag, inside, z, coords, q_all, opts, u, error)
      !! The field the atoms marked `inside` sit in, over `mol`'s basis
      !!
      !! Works for a monomer and a dimer alike -- `inside` and `near` are what
      !! change. Comes back unallocated when there is no field, which is how the
      !! caller knows the SCF energy needs no correction.
      !!
      !! Two kinds of outside atom contribute, and the split is what `resppc`
      !! decides. A near fragment gives its bare nuclei here and its electrons
      !! through the exact Coulomb operator below. A distant one gives a single
      !! partial charge per atom, which already carries both -- that is what
      !! `q_A = Z_A - population_A` means, and why the point-charge form is an
      !! approximation to the electron term alone rather than to the whole
      !! interaction.
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: group_z(:)
      character(len=2), intent(in) :: group_sym(:)
      real(dp), intent(in) :: group_xyz(:, :)
      integer, intent(in) :: near(:)
      type(fragment_t), intent(in) :: frag(:)
      integer, intent(in) :: n_frag
      logical, intent(in) :: inside(:)
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords(:, :)
      real(dp), allocatable, intent(in) :: q_all(:)
      type(fmo_options_t), intent(in) :: opts
      real(dp), allocatable, intent(out) :: u(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: matrices(:, :, :), points(:, :), weight(:), j_near(:, :)
      logical, allocatable :: is_near(:)
      integer :: n_atoms, i, k, g, n_out

      if (opts%esp == "none") return
      n_atoms = size(z)
      n_out = count(.not. inside)
      if (n_out == 0) return

      ! Which outside atoms belong to a fragment being treated exactly.
      allocate (is_near(n_atoms), source=.false.)
      do k = 1, size(near)
         is_near(frag(near(k))%atoms) = .true.
      end do

      ! Distant atoms drop out entirely when they are being ignored, so the
      ! integrals for them are never evaluated rather than evaluated and scaled
      ! by nothing.
      allocate (points(3, n_out), weight(n_out))
      g = 0
      do i = 1, n_atoms
         if (inside(i)) cycle
         if (.not. is_near(i) .and. opts%far_field == "ignore") cycle
         g = g + 1
         points(:, g) = coords(:, i)
         if (is_near(i)) then
            weight(g) = real(z(i), dp)      !! nucleus only; electrons via J below
         else
            weight(g) = q_all(i)            !! nucleus and electrons, approximated
         end if
      end do

      ! Nothing near and nothing distant that counts: there is no field, and
      ! saying so with an unallocated result spares the caller a zero matrix it
      ! would have to add to everything anyway.
      if (g == 0 .and. size(near) == 0) return

      allocate (u(mol%nao, mol%nao), source=0.0_dp)
      if (g > 0) then
         ! An electron carries charge -1, so a positive charge lowers its
         ! energy: the operator is -sum_g w_g/|r - R_g|.
         call esp_matrices(mol, points(:, 1:g), matrices, error)
         if (error%has_error()) return
         do i = 1, g
            u = u - weight(i)*matrices(:, :, i)
         end do
      end if
      if (size(near) == 0) return

      call local_coulomb(frag, group_z, group_sym, group_xyz, mol%nao, near, opts, &
                         j_near, error)
      if (error%has_error()) return
      u = u + j_near
   end subroutine embedding_operator

   subroutine pair_term(frag, n_frag, a, b, z, coords, q_all, opts, res, &
                        all_converged, error)
      !! One dimer, in the field of every fragment outside it
      type(fragment_t), intent(in) :: frag(:)
      integer, intent(in) :: n_frag, a, b
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords(:, :)
      real(dp), allocatable, intent(in) :: q_all(:)
      type(fmo_options_t), intent(in) :: opts
      type(fmo_result_t), intent(inout) :: res
      logical, intent(inout) :: all_converged
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: bounds(:, :), d_split(:, :), u(:, :)
      real(dp), allocatable :: pair_xyz(:, :)
      logical, allocatable :: inside(:)
      integer, allocatable :: near(:)
      integer, allocatable :: pair_z(:)
      character(len=2), allocatable :: pair_sym(:)
      real(dp) :: e_internal, e_resp, e_pair
      integer :: na, nb

      na = frag(a)%nao
      nb = frag(b)%nao

      pair_z = [frag(a)%z, frag(b)%z]
      pair_sym = [frag(a)%sym, frag(b)%sym]
      pair_xyz = reshape([frag(a)%xyz, frag(b)%xyz], [3, size(pair_z)])
      call build_libcint_molecule(pair_z, pair_sym, pair_xyz, trim(opts%basis), mol, error)
      if (error%has_error()) return

      ! The dimer is fragment a's atoms then fragment b's, and libcint orders
      ! functions by atom, so the monomer blocks sit contiguously front and
      ! back. Worth checking rather than assuming: a silent mismatch here would
      ! be a plausible-looking wrong answer rather than a failure.
      if (mol%nao /= na + nb) then
         call error%set(ERROR_VALIDATION, "fmo: the dimer basis is not the two monomer "// &
                        "bases end to end ("//to_char(na)//" + "//to_char(nb)//" /= "// &
                        to_char(mol%nao)//")")
         return
      end if
      call schwarz_bounds(mol, bounds, error)
      if (error%has_error()) return

      ! The two monomer densities side by side: what the dimer's own Coulomb
      ! contribution is subtracted with, and what dD is measured against.
      allocate (d_split(mol%nao, mol%nao), source=0.0_dp)
      d_split(1:na, 1:na) = frag(a)%density
      d_split(na + 1:, na + 1:) = frag(b)%density

      allocate (inside(size(z)), source=.false.)
      inside(frag(a)%atoms) = .true.
      inside(frag(b)%atoms) = .true.

      ! The pair's own near set, not the union of the two monomers': a fragment
      ! near either one is near the pair.
      call near_fragments(frag, n_frag, [a, b], effective_resppc(opts), near, error)
      if (error%has_error()) return

      call embedding_operator(mol, pair_z, pair_sym, pair_xyz, near, frag, n_frag, &
                              inside, z, coords, q_all, opts, u, error)
      if (error%has_error()) return

      if (allocated(u)) then
         call run_libcint_rhf(mol, frag(a)%nelec + frag(b)%nelec, opts%scf_max_iter, &
                              opts%scf_energy_tol, opts%scf_density_tol, .false., &
                              scf, error, h_extra=u)
      else
         call run_libcint_rhf(mol, frag(a)%nelec + frag(b)%nelec, opts%scf_max_iter, &
                              opts%scf_energy_tol, opts%scf_density_tol, .false., &
                              scf, error)
      end if
      if (error%has_error()) return
      if (.not. scf%converged) all_converged = .false.

      e_internal = scf%energy
      e_resp = 0.0_dp
      if (allocated(u)) then
         e_internal = e_internal - sum(scf%density*u)
         e_resp = sum((scf%density - d_split)*u)
      end if

      if (opts%expansion == "mbe") then
         ! Total embedded energies, straight difference, nothing else.
         e_pair = scf%energy - frag(a)%energy_total - frag(b)%energy_total
         e_resp = 0.0_dp
      else
         e_pair = e_internal - frag(a)%energy - frag(b)%energy
      end if

      res%pair_correction(a, b) = e_pair + e_resp
      res%pair_sum = res%pair_sum + e_pair
      res%response_sum = res%response_sum + e_resp
   end subroutine pair_term

   subroutine near_fragments(frag, n_frag, group, resppc, near, error)
      !! Which fragments outside `group` are close enough to need the exact term
      !!
      !! FMO's separation measure: the smallest interatomic distance between the
      !! two fragments, over the sum of those two atoms' van der Waals radii.
      !! Contact distance, not centre-to-centre, so an extended fragment counts
      !! as near if any part of it is.
      type(fragment_t), intent(in) :: frag(:)
      integer, intent(in) :: n_frag
      integer, intent(in) :: group(:)
      real(dp), intent(in) :: resppc
      integer, allocatable, intent(out) :: near(:)
      type(error_t), intent(inout) :: error

      logical, allocatable :: keep(:)
      real(dp) :: r, best
      integer :: k, g, ia, ib

      allocate (keep(n_frag), source=.false.)
      do k = 1, n_frag
         if (any(group == k)) cycle
         if (resppc < 0.0_dp) then
            keep(k) = .true.          !! approximation disabled: everything exact
            cycle
         end if

         best = huge(1.0_dp)
         do g = 1, size(group)
            do ia = 1, size(frag(group(g))%atoms)
               do ib = 1, size(frag(k)%atoms)
                  call unitless_distance(frag(group(g)), ia, frag(k), ib, r, error)
                  if (error%has_error()) return
                  best = min(best, r)
               end do
            end do
         end do
         keep(k) = best <= resppc
      end do

      near = pack([(k, k=1, n_frag)], keep)
   end subroutine near_fragments

   subroutine unitless_distance(fa, ia, fb, ib, r, error)
      !! |R_A - R_B| / (vdw_A + vdw_B), the separation FMO measures in
      real(dp), parameter :: ANGSTROM_TO_BOHR = 1.8897261254578281_dp
      type(fragment_t), intent(in) :: fa, fb
      integer, intent(in) :: ia, ib
      real(dp), intent(out) :: r
      type(error_t), intent(inout) :: error

      real(dp) :: scale

      scale = (element_vdw_radius(fa%z(ia)) + element_vdw_radius(fb%z(ib)))*ANGSTROM_TO_BOHR
      if (scale <= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "fmo: no van der Waals radius for element "// &
                        to_char(fa%z(ia))//" or "//to_char(fb%z(ib))//", so the "// &
                        "separation the point-charge cutoff is measured in is undefined")
         return
      end if
      r = norm2(fa%xyz(:, ia) - fb%xyz(:, ib))/scale
   end subroutine unitless_distance

   subroutine local_coulomb(frag, group_z, group_sym, group_xyz, group_nao, near, &
                            opts, j_near, error)
      !! The exact Coulomb operator of the near fragments, on the group's basis
      !!
      !! Built over a molecule holding the group and its near neighbours only,
      !! with the group's own density zeroed -- so what comes back on the
      !! group's block is `sum_{K near} sum_{ls in K} D^K_ls (mn|ls)` and
      !! nothing else. No subtraction, and no dependence on the size of the
      !! system beyond the neighbourhood, which is the point of having a cutoff
      !! at all.
      type(fragment_t), intent(in) :: frag(:)
      integer, intent(in) :: group_z(:)
      character(len=2), intent(in) :: group_sym(:)
      real(dp), intent(in) :: group_xyz(:, :)
      integer, intent(in) :: group_nao
      integer, intent(in) :: near(:)
      type(fmo_options_t), intent(in) :: opts
      real(dp), allocatable, intent(out) :: j_near(:, :)
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: local
      type(direct_stats_t) :: stats
      integer, allocatable :: z(:)
      character(len=2), allocatable :: sym(:)
      real(dp), allocatable :: xyz(:, :), d(:, :), zero_h(:, :), bounds(:, :), j_full(:, :)
      integer :: k, at, nao_k, expect

      allocate (j_near(group_nao, group_nao), source=0.0_dp)
      if (size(near) == 0) return

      z = group_z
      sym = group_sym
      xyz = group_xyz
      expect = group_nao
      do k = 1, size(near)
         z = [z, frag(near(k))%z]
         sym = [sym, frag(near(k))%sym]
         xyz = reshape([xyz, frag(near(k))%xyz], [3, size(z)])
         expect = expect + frag(near(k))%nao
      end do

      call build_libcint_molecule(z, sym, xyz, trim(opts%basis), local, error)
      if (error%has_error()) return
      if (local%nao /= expect) then
         call error%set(ERROR_VALIDATION, "fmo: the local supersystem has "// &
                        to_char(local%nao)//" basis functions where its parts have "// &
                        to_char(expect)//"; the atom ordering is not what the block "// &
                        "layout assumes")
         return
      end if
      call schwarz_bounds(local, bounds, error)
      if (error%has_error()) return

      ! Only the neighbours carry density. The group's own block stays zero, so
      ! its own Coulomb contribution never enters and needs no removing.
      allocate (d(local%nao, local%nao), source=0.0_dp)
      at = group_nao
      do k = 1, size(near)
         nao_k = frag(near(k))%nao
         d(at + 1:at + nao_k, at + 1:at + nao_k) = frag(near(k))%density
         at = at + nao_k
      end do

      allocate (zero_h(local%nao, local%nao), source=0.0_dp)
      allocate (j_full(local%nao, local%nao))
      call build_fock_direct(local, zero_h, d, bounds, j_full, stats, error, &
                             k_scale=0.0_dp, j_scale=1.0_dp)
      if (error%has_error()) return
      j_near = j_full(1:group_nao, 1:group_nao)
   end subroutine local_coulomb

   subroutine solve_fragment(frag, n_frag, which, z, coords, q_all, opts, &
                             all_converged, error, bare)
      !! One monomer: build it, field it, solve it, read its charges, drop it
      !!
      !! The whole of a fragment's work for one outer pass, and the unit a rank
      !! would be handed. Its molecule exists only inside this call, which is
      !! what lets the caller hold nothing heavy and lets two ranks do two
      !! fragments without sharing anything but the geometry they both started
      !! with.
      type(fragment_t), intent(inout) :: frag(:)
      integer, intent(in) :: n_frag, which
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords(:, :)
      real(dp), allocatable, intent(in) :: q_all(:)
      type(fmo_options_t), intent(in) :: opts
      logical, intent(inout) :: all_converged
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: bare
         !! Solve it in vacuum. True only for the very first pass, which has no
         !! field to be solved against because no fragment has a density yet --
         !! the field is what the passes after it are for.

      type(libcint_molecule_t) :: mol
      real(dp), allocatable :: bounds(:, :), u(:, :), q(:)
      logical, allocatable :: inside(:)
      logical :: isolated

      isolated = .false.
      if (present(bare)) isolated = bare

      call open_fragment(frag(which)%z, frag(which)%sym, frag(which)%xyz, opts, &
                         mol, bounds, error)
      if (error%has_error()) return

      if (.not. isolated) then
         allocate (inside(size(z)), source=.false.)
         inside(frag(which)%atoms) = .true.
         call embedding_operator(mol, frag(which)%z, frag(which)%sym, frag(which)%xyz, &
                                 frag(which)%near, frag, n_frag, inside, z, coords, &
                                 q_all, opts, u, error)
         if (error%has_error()) return
      end if

      call inner_scf(frag(which), mol, opts, error, u, all_converged)
      if (error%has_error()) return

      if (opts%esp /= "none") then
         call fragment_charges(mol, frag(which)%density, opts%far_field, q, error)
         if (error%has_error()) return
         frag(which)%charges = q
      end if
   end subroutine solve_fragment

   subroutine calculate_monomers(frag, n_frag, z, coords, opts, res, all_converged, error, comm)
      !! The outer SCF: every monomer in the field of all the others, iterated
      !!
      !! One pass solves every fragment against the field the previous pass
      !! left. Converged when the monomer energies stop moving, which is when
      !! the densities and the field they make agree.
      !!
      !! Every fragment within a pass is independent -- they all read the
      !! previous pass's densities, none reads this pass's -- so a pass is a
      !! bag of independent tasks with a barrier after it. That is the shape a
      !! distributed version needs: hand the pass out, wait, share the
      !! densities, decide whether to go again.
      type(fragment_t), intent(inout) :: frag(:)
      integer, intent(in) :: n_frag
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords(:, :)
      type(fmo_options_t), intent(in) :: opts
      type(fmo_result_t), intent(inout) :: res
      logical, intent(inout) :: all_converged
      type(error_t), intent(inout) :: error
      type(comm_t), intent(in), optional :: comm

      real(dp), allocatable :: q_all(:)
      real(dp) :: e_sum, e_prev
      integer :: i, outer

      ! Isolated fragments, to start the field from something.
      call all_charges(frag, n_frag, size(z), opts, q_all, error)
      if (error%has_error()) return
      do i = 1, n_frag
         if (.not. mine(i, comm)) cycle
         call solve_fragment(frag, n_frag, i, z, coords, q_all, opts, all_converged, &
                             error, bare=.true.)
         if (error%has_error()) return
      end do
      call exchange_monomers(frag, n_frag, comm)

      if (opts%esp == "none") then
         res%converged = .true.
         res%outer_iterations = 1
         return
      end if

      e_prev = sum(frag(:)%energy)
      do outer = 1, opts%max_outer
         call all_charges(frag, n_frag, size(z), opts, q_all, error)
         if (error%has_error()) return

         ! Independent within a pass: every fragment reads the densities the
         ! last pass left and none reads this pass's, so who computes which is
         ! free. The exchange after is the barrier.
         do i = 1, n_frag
            if (.not. mine(i, comm)) cycle
            call solve_fragment(frag, n_frag, i, z, coords, q_all, opts, &
                                all_converged, error)
            if (error%has_error()) return
         end do
         call exchange_monomers(frag, n_frag, comm)

         e_sum = sum(frag(:)%energy)
         res%outer_iterations = outer
         res%outer_change = abs(e_sum - e_prev)
         call logger%verbose("  fmo outer "//to_char(outer)//": monomer sum "// &
                             to_char(e_sum)//", moved "//to_char(res%outer_change))
         if (res%outer_change < opts%outer_tol) then
            res%converged = .true.
            return
         end if
         e_prev = e_sum
      end do

      call error%set(ERROR_VALIDATION, "fmo: the outer SCF did not settle in "// &
                     to_char(opts%max_outer)//" passes; the monomer sum was still "// &
                     "moving by "//to_char(res%outer_change)//" Hartree")
   end subroutine calculate_monomers

   subroutine calculate_polymers(frag, n_frag, z, coords, opts, res, all_converged, error, comm)
      !! Every pair, in the field the converged monomers make
      !!
      !! Independent of each other and of everything else once the monomers have
      !! settled, so this is a single bag of tasks with no barrier inside it --
      !! the easy half to distribute, and the expensive one, being quadratic in
      !! the fragment count where the monomer loop is linear.
      type(fragment_t), intent(inout) :: frag(:)
      integer, intent(in) :: n_frag
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords(:, :)
      type(fmo_options_t), intent(in) :: opts
      type(fmo_result_t), intent(inout) :: res
      logical, intent(inout) :: all_converged
      type(error_t), intent(inout) :: error
      type(comm_t), intent(in), optional :: comm

      real(dp), allocatable :: q_all(:), totals(:)
      integer :: i, j, pair_index

      call all_charges(frag, n_frag, size(z), opts, q_all, error)
      if (error%has_error()) return

      ! One bag of tasks, no barrier inside it: nothing here feeds anything
      ! else here. Quadratic in the fragment count where the monomer loop is
      ! linear, so this is the half worth spreading.
      pair_index = 0
      do i = 1, n_frag
         do j = i + 1, n_frag
            pair_index = pair_index + 1
            if (.not. mine(pair_index, comm)) cycle
            call pair_term(frag, n_frag, i, j, z, coords, q_all, opts, res, &
                           all_converged, error)
            if (error%has_error()) return
         end do
      end do

      if (spread_over(comm)) then
         ! Each rank accumulated only its own pairs, so the totals are partial
         ! sums and add up to the whole.
         allocate (totals(2 + n_frag*n_frag))
         totals(1) = res%pair_sum
         totals(2) = res%response_sum
         totals(3:) = reshape(res%pair_correction, [n_frag*n_frag])
         call allreduce(comm, totals, size(totals), MPI_SUM)
         res%pair_sum = totals(1)
         res%response_sum = totals(2)
         res%pair_correction = reshape(totals(3:), [n_frag, n_frag])
      end if
   end subroutine calculate_polymers

   function spread_over(comm) result(many)
      !! Whether there is more than one rank to spread over
      type(comm_t), intent(in), optional :: comm
      logical :: many

      many = .false.
      if (.not. present(comm)) return
      many = comm%size() > 1
   end function spread_over

   function mine(task, comm) result(owned)
      !! Whether this rank owns a task, round robin
      !!
      !! Static rather than handed out on demand, because FMO's tasks are all
      !! much of a size -- every monomer is one fragment, every pair is two --
      !! so the imbalance a task server exists to absorb is not there to absorb.
      integer, intent(in) :: task
      type(comm_t), intent(in), optional :: comm
      logical :: owned

      owned = .true.
      if (.not. spread_over(comm)) return
      owned = mod(task - 1, comm%size()) == comm%rank()
   end function mine

   subroutine exchange_monomers(frag, n_frag, comm)
      !! Share what each rank computed this pass with every other
      !!
      !! A sum-reduce over buffers that are zero where a rank computed nothing,
      !! which makes it a gather without needing displacements for fragments of
      !! different sizes. Densities, energies and charges only -- small, and the
      !! only things the next pass reads.
      type(fragment_t), intent(inout) :: frag(:)
      integer, intent(in) :: n_frag
      type(comm_t), intent(in), optional :: comm

      real(dp), allocatable :: buf(:)
      integer :: f, at, n, total

      if (.not. spread_over(comm)) return

      total = 0
      do f = 1, n_frag
         total = total + frag(f)%nao*frag(f)%nao + 2 + size(frag(f)%atoms)
      end do
      allocate (buf(total), source=0.0_dp)

      at = 0
      do f = 1, n_frag
         n = frag(f)%nao*frag(f)%nao
         if (mine(f, comm)) then
            if (allocated(frag(f)%density)) buf(at + 1:at + n) = reshape(frag(f)%density, [n])
            buf(at + n + 1) = frag(f)%energy
            buf(at + n + 2) = frag(f)%energy_total
            if (allocated(frag(f)%charges)) then
               buf(at + n + 3:at + n + 2 + size(frag(f)%atoms)) = frag(f)%charges
            end if
         end if
         at = at + n + 2 + size(frag(f)%atoms)
      end do

      call allreduce(comm, buf, size(buf), MPI_SUM)

      at = 0
      do f = 1, n_frag
         n = frag(f)%nao*frag(f)%nao
         if (.not. allocated(frag(f)%density)) allocate (frag(f)%density(frag(f)%nao, frag(f)%nao))
         frag(f)%density = reshape(buf(at + 1:at + n), [frag(f)%nao, frag(f)%nao])
         frag(f)%energy = buf(at + n + 1)
         frag(f)%energy_total = buf(at + n + 2)
         if (.not. allocated(frag(f)%charges)) allocate (frag(f)%charges(size(frag(f)%atoms)))
         frag(f)%charges = buf(at + n + 3:at + n + 2 + size(frag(f)%atoms))
         at = at + n + 2 + size(frag(f)%atoms)
      end do
   end subroutine exchange_monomers

   subroutine inner_scf(f, mol, opts, error, u, all_converged)
      !! The inner SCF: this fragment's orbitals, against a fixed external field
      type(fragment_t), intent(inout) :: f
      type(libcint_molecule_t), intent(in) :: mol
      type(fmo_options_t), intent(in) :: opts
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(in), optional :: u(:, :)
      logical, intent(inout), optional :: all_converged

      type(rhf_result_t) :: scf
      logical :: embedded

      embedded = .false.
      if (present(u)) embedded = allocated(u)

      if (embedded) then
         call run_libcint_rhf(mol, f%nelec, opts%scf_max_iter, opts%scf_energy_tol, &
                              opts%scf_density_tol, .false., scf, error, h_extra=u)
      else
         call run_libcint_rhf(mol, f%nelec, opts%scf_max_iter, opts%scf_energy_tol, &
                              opts%scf_density_tol, .false., scf, error)
      end if
      if (error%has_error()) return
      if (present(all_converged)) then
         if (.not. scf%converged) all_converged = .false.
      end if

      ! The internal energy: what the SCF reported, less its interaction with
      ! the field. `h_extra` enters H linearly, so that interaction is exactly
      ! Tr(D u) and nothing else has to be unpicked.
      f%energy_total = scf%energy
      f%energy = scf%energy
      if (embedded) f%energy = f%energy - sum(scf%density*u)
      f%density = scf%density
   end subroutine inner_scf

   subroutine fragment_charges(mol, density, scheme, q, error)
      !! Atomic charges for a fragment whose density is already converged
      !!
      !! Taken while the molecule is in hand rather than rebuilt for it later:
      !! the overlap and the ESP grid both need it, and it is about to be thrown
      !! away.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      character(len=*), intent(in) :: scheme
      real(dp), allocatable, intent(out) :: q(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: s(:, :)

      if (trim(scheme) == "chelpg") then
         call chelpg_charges(mol, density, q, error)
      else
         call mol%overlap(s)
         call mulliken_charges(mol, density, s, q, error)
      end if
   end subroutine fragment_charges

   subroutine report_charges(frag, n_frag, n_atoms, charges, error)
      !! Mulliken charges over the whole system, for the result record
      type(fragment_t), intent(in) :: frag(:)
      integer, intent(in) :: n_frag, n_atoms
      real(dp), allocatable, intent(out) :: charges(:)
      type(error_t), intent(inout) :: error

      integer :: f

      allocate (charges(n_atoms), source=0.0_dp)
      do f = 1, n_frag
         if (.not. allocated(frag(f)%charges)) cycle
         charges(frag(f)%atoms) = frag(f)%charges
      end do
   end subroutine report_charges

end module mqc_libcint_fmo
