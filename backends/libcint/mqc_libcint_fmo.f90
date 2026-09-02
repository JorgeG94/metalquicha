!! The fragment molecular orbital method, to any order
module mqc_libcint_fmo
   !! FMO_n: two nested self-consistencies, then a pass over every n-mer.
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
   !! with its polarised density, not counting its interaction with the field --
   !! every group of fragments `S` contributes what its own subsets did not:
   !!
   !!     dE_S = E'_S + Tr(dD_S u_S) - sum over proper subsets T of S of dE_T
   !!
   !!     E = sum over all S up to the level of dE_S
   !!
   !! For a pair that unrolls to `E'_IJ - E'_I - E'_J` plus the pair's density
   !! response to the field outside it, and for a trimer to the usual three-body
   !! expression, neither being written down anywhere. `dD_S` is the n-mer's
   !! density less its members' densities laid side by side.
   !!
   !! The response belongs inside the recursion. Left outside it, a pair's
   !! response is never cancelled by a trimer containing that pair, and the
   !! total misses exactness by precisely the sum of them -- invisible at level
   !! two, where nothing contains a pair.
   !!
   !! **Level `n` on `n` fragments is not an approximation.** The corrections
   !! telescope to the supermolecular energy, the same way FMO2 on two fragments
   !! does. `check_fmo` climbs the ladder and ends on that equality rather than
   !! on a tolerance, which is what pins the many-body coefficients down: a
   !! wrong one still looks plausible at level two and cannot survive being
   !! asked to cancel exactly.
   !!
   !! **Whole molecules by default, and a covalent cut has to be asked for.**
   !! `bond_breaking = "none"` refuses a partition that severs a bond rather
   !! than answering it; `"afo"` detaches the bond with an adjusted frozen
   !! orbital instead, and is restricted to `esp = "none"` for now.
   !!
   !! The refusal is not a formality, which is why it stayed the default.
   !! Cutting a single bond leaves both fragments with an odd electron count,
   !! which the closed-shell check would catch on its own; but cut an even
   !! number per fragment -- a ring, a double bond -- and every count stays
   !! even. Cyclopropane split into three CH2 used to come back 0.28 Hartree
   !! low, which is 176 kcal/mol in the shape of an answer. Connectivity is
   !! therefore checked directly, with the criterion [[mqc_bond_perception]]
   !! uses everywhere else.
   !!
   !! When a bond *is* detached, the orbital that stands in for it comes from
   !! [[mqc_libcint_afo]] and is held by [[mqc_fock_projector]]. What belongs
   !! here is only the assembly: `assemble_group` decides a group's boundaries
   !! from its own members, every time, because a bond cut between two monomers
   !! is whole again inside the dimer holding both ends. Inheriting that
   !! decision from the members instead is what made an earlier capped version
   !! 11 Hartree wrong.
   !!
   !! **Cost.** There are C(N,n) n-mers, so level three on twenty fragments is
   !! 1140 SCFs against 190 for level two. Nothing here refuses a high level --
   !! the expansion is generic and a cap would be arbitrary -- but the binomial
   !! is the whole story and it is not gentle.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger, verbose_level, debug_level, info_level
   use pic_io, only: to_char
   use mqc_convergence_report, only: convergence_header, convergence_footer
   use pic_mpi_lib, only: comm_t, allreduce, MPI_SUM
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_elements, only: element_vdw_radius
   use mqc_physical_constants, only: ANGSTROM_TO_BOHR
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_bond_perception, only: connected_components, find_severed_bonds, severed_bond_t
   use mqc_libcint_afo, only: afo_model_t, afo_options_t, afo_hybrid_t, build_afo_model, &
                              bond_hybrid, cuts_outside_group, group_electron_shift, &
                              build_group_frozen
   use mqc_fock_projector, only: fock_projector_t, build_frozen_basis
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule, &
                                    atom_ao_blocks
   use mqc_libcint_direct, only: schwarz_bounds, build_fock_direct, direct_stats_t
   use mqc_libcint_esp, only: esp_matrices
   use mqc_libcint_charges, only: mulliken_charges, chelpg_charges
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_scf_types, only: scf_numerics_t
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
      character(len=16) :: bond_breaking = "none"
         !! How a cut covalent bond is represented.
         !!
         !! `"none"` refuses a partition that cuts one, naming the two atoms
         !! and the two fragments they were put in. That is the default and
         !! it is not a formality: cutting an even number of bonds per
         !! fragment leaves every electron count even, so nothing else
         !! objects, and cyclopropane split into three CH2 used to come back
         !! 0.28 Hartree low.
         !!
         !! `"afo"` detaches the bond with an adjusted frozen orbital. A model
         !! system around the bond is solved and localized, the orbital on the
         !! bond is reduced to the detached atom's own functions, and that
         !! hybrid is then frozen -- empty in the fragment that gets nothing of
         !! the bond, occupied in the one that gets all of it. See
         !! [[mqc_libcint_afo]].
         !!
         !! Only sound without an embedding field so far, so `"afo"` requires
         !! `esp = "none"` and says so. A frozen orbital and a field both
         !! describe the bond region -- the field already supplies the
         !! neighbour's nucleus and density where the orbital supplies the
         !! bond -- so the detached atom's share has to come out of the field
         !! first. That is clean for point charges and not defined for an
         !! exact density.
      real(dp) :: cap_scale = 1.0_dp
         !! Where a hydrogen cap sits along the bond it closes, for the
         !! many-body path; see [[mqc_physical_fragment]]. Not used by
         !! `"afo"`, which caps only the model systems it builds and derives
         !! their scale from covalent radii rather than from a deck.
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
      integer :: level = 2
         !! How many fragments at a time. Two is FMO2, three is FMO3, and the
         !! expansion is truncated there.
         !!
         !! The cost is the binomial: with N fragments there are C(N,n) n-mers,
         !! so level 3 on twenty fragments is 1140 SCFs against 190 for level 2,
         !! and level 7 is not a good idea on anything. It is allowed because
         !! the expansion is generic and refusing it would be arbitrary, not
         !! because it is advisable.
         !!
         !! Level equal to the fragment count is the whole expansion and so is
         !! exact -- the corrections telescope to the supermolecular energy.
         !! That is a useful thing to be able to ask for, and `check_fmo`
         !! asserts it at every level a small system allows.
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
      type(scf_numerics_t) :: scf
         !! How each fragment SCF is driven. The three fields below stay
         !! per-fragment -- a fragment is a smaller problem and gets its own
         !! budget and tolerances -- but the accelerator, DIIS subspace, level
         !! shift and linear-dependence threshold are properties of the
         !! calculation, and none of them reached a fragment before.
         !! **Only the drive settings here are read** -- the accelerator, DIIS
         !! subspace, level shift, linear-dependence threshold and incremental
         !! Fock switch. Its `max_iter`, `energy_tol` and `density_tol` are
         !! NOT: the three bare fields below are, and they are passed
         !! positionally so they win. Two declarations of one concept in one
         !! type is a trap, and this comment is the guard rail until the bare
         !! three are folded in -- which cannot happen until their deliberately
         !! tighter values have somewhere else to live.
      integer :: scf_max_iter = 100
      real(dp) :: scf_energy_tol = 1.0e-9_dp
      real(dp) :: scf_density_tol = 1.0e-7_dp
   end type fmo_options_t

   type :: fmo_result_t
      !! The energy, and enough of the parts to see where it came from
      real(dp) :: energy = 0.0_dp                !! The FMO2 total
      real(dp) :: monomer_sum = 0.0_dp           !! sum_I E'_I
      real(dp) :: pair_sum = 0.0_dp              !! sum of every n-mer correction
      real(dp) :: response_sum = 0.0_dp          !! sum of Tr(dD u), the last term
      integer :: outer_iterations = 0            !! passes of the monomer SCF
      real(dp) :: outer_change = 0.0_dp          !! last movement of the monomer sum
      logical :: converged = .false.
      real(dp), allocatable :: monomer_energy(:)     !! E'_I
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
      integer :: n_caps = 0
         !! Hydrogen caps closing cut bonds, held at the end of `z`, `sym`
         !! and `xyz` and deliberately absent from `atoms`. `atoms` maps a
         !! fragment's entries back onto the system and a cap answers to no
         !! system atom, so it must not appear there -- it is used as a
         !! scatter index in several places and a zero would write outside
         !! the array. The consequence is that `size(atoms)` counts real
         !! atoms and `size(z)` counts what the SCF actually sees.
      integer, allocatable :: near(:)            !! fragments close enough for the exact term
   end type fragment_t

   !> Everything the system's cut bonds imply, worked out once
   type :: afo_context_t
      !! Built before any fragment is solved, because a hybrid is a property of
      !! the bond and its surroundings and not of whoever is being solved. One
      !! model system per cut bond, whatever the fragment count.
      logical :: active = .false.
      integer :: n_cuts = 0
      type(severed_bond_t), allocatable :: cuts(:)
      type(afo_hybrid_t), allocatable :: hybrid(:)
      character(len=2), allocatable :: sym(:)
         !! System-wide symbols. A group holding the attached end of a bond has
         !! to name the detached atom to ghost it, and that atom is one it does
         !! not own, so its symbol cannot come from any fragment.
   end type afo_context_t

   !> One n-mer as it will be handed to an SCF, boundaries included
   type :: group_t
      integer :: n_real = 0    !! Atoms the group owns; ghosts follow them
      integer :: nelec = 0     !! With the boundary shift already applied
      integer, allocatable :: z(:)
      character(len=2), allocatable :: sym(:)
      real(dp), allocatable :: xyz(:, :)
      logical, allocatable :: ghost(:)
      integer :: n_bound = 0   !! Cut bonds this group is still cut across
      integer, allocatable :: bda_slot(:)  !! Where each boundary's BDA sits here
      logical, allocatable :: occupied(:)  !! True when this group gets the bond
      integer, allocatable :: cut_of(:)    !! Which system cut each boundary is
   end type group_t

   !> Where a frozen virtual is held, in Hartree.
   !>
   !> Modest on purpose. The blocks are decoupled rather than penalised, so this
   !> has only to lift the frozen virtuals clear of the occupied manifold, and a
   !> larger one costs precision -- the back transform spreads it over every
   !> element, leaving the unfrozen block clean only to about `shift * epsilon`.
   real(dp), parameter :: AFO_SHIFT = 1.0e3_dp

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
      type(afo_context_t) :: afo
      integer :: n_atoms, n_frag, i
      logical :: all_converged

      n_atoms = size(atomic_numbers)
      if (size(owner) /= n_atoms .or. size(coordinates, 2) /= n_atoms) then
         call error%set(ERROR_VALIDATION, "fmo: owner and coordinates must cover "// &
                        "every atom")
         return
      end if

      call build_fragments(atomic_numbers, symbols, coordinates, owner, opts, &
                           frag, n_frag, afo, error, comm)
      if (error%has_error()) return

      allocate (res%monomer_energy(n_frag), source=0.0_dp)
      all_converged = .true.

      call calculate_monomers(frag, n_frag, atomic_numbers, coordinates, opts, afo, res, &
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

      call calculate_polymers(frag, n_frag, atomic_numbers, coordinates, opts, afo, res, &
                              all_converged, error, comm)
      if (error%has_error()) return

      ! `response_sum` is inside `pair_sum` already; it is reported, not added.
      res%energy = res%monomer_sum + res%pair_sum

      ! Reported, not used. What the fragments look like once the field has
      ! settled is the first thing to look at when a number seems wrong.
      call report_charges(frag, n_frag, n_atoms, res%charges, error)
      if (error%has_error()) return

      res%converged = res%converged .and. all_converged
      if (.not. all_converged) then
         ! Refused unless the deck said otherwise. `allow_crap_scf` means the
         ! same here as everywhere else -- keep a total built partly from
         ! unconverged pieces, as a result to follow up rather than to trust --
         ! and it was not reaching this path at all, so a fragmented job that a
         ! deck had explicitly allowed to finish was refused anyway.
         if (.not. opts%scf%allow_crap_scf) then
            call error%set(ERROR_VALIDATION, "fmo: at least one fragment SCF did not "// &
                           "converge, so the total is not trustworthy. Set "// &
                           "keywords.scf.allow_crap_scf to finish anyway.")
            return
         end if
         call logger%warning("  fmo: at least one fragment SCF did not converge, and "// &
                             "allow_crap_scf kept the run. The total is built partly "// &
                             "from unconverged fragments.")
      end if
   end subroutine run_fmo2

   subroutine fmo_guess_kind(opts, kind)
      !! The deck's initial guess for a fragment SCF, as an `SCF_GUESS_*` kind
      !!
      !! Fragment SCFs never received one: `keywords.scf.guess` reached the
      !! expansion and stopped there, so every fragment used whatever the SCF
      !! resolves for itself no matter what the deck asked for.
      !!
      !! Absent rather than a value when the deck said nothing. `auto` means
      !! "the backend picks", and forwarding a resolved kind for it would pin
      !! the choice here instead -- a behaviour change for every deck that
      !! never mentioned a guess. An unparseable spelling falls back the same
      !! way: a guess is a starting point, and refusing a whole fragmented
      !! calculation over one is out of proportion.
      !!
      !! A subroutine filling a local, not a function: an unallocated
      !! allocatable VARIABLE arrives absent at an optional dummy, but the same
      !! thing as a function RESULT segfaults -- found the hard way, on
      !! check_fmo_mpi.
      use mqc_libcint_atomic_guess, only: parse_guess_name
      use mqc_libcint_rhf, only: SCF_GUESS_SAD, SCF_GUESS_SAC, SCF_GUESS_PROJ
      type(fmo_options_t), intent(in) :: opts
      integer, allocatable, intent(out) :: kind

      integer :: parsed
      type(error_t) :: guess_error

      if (trim(opts%scf%guess) == "auto" .or. len_trim(opts%scf%guess) == 0) return
      call parse_guess_name(opts%scf%guess, parsed, guess_error)
      if (guess_error%has_error()) then
         call guess_error%clear()
         return
      end if

      ! `sad`, `sac` and `projection` need a `guess_density` alongside the kind,
      ! and a fragment SCF has nobody to build one for it -- `run_libcint_rhf`
      ! refuses outright when the kind arrives without the density. Forwarding
      ! them would turn the most natural thing a user writes into a fatal error
      ! on every fragment, which is the opposite of what this routine is for, so
      ! they fall back to the backend's own choice and say so once.
      if (parsed == SCF_GUESS_SAD .or. parsed == SCF_GUESS_SAC .or. &
          parsed == SCF_GUESS_PROJ) then
         call logger%warning("  the '"//trim(opts%scf%guess)//"' initial guess needs a "// &
                             "guess density, which a fragment SCF cannot be given yet; "// &
                             "fragments use the backend's default guess instead")
         return
      end if
      allocate (kind, source=parsed)
   end subroutine fmo_guess_kind

   subroutine open_fragment(frag_z, frag_sym, frag_xyz, opts, mol, bounds, error, ghost)
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
      logical, intent(in), optional :: ghost(:)
         !! Atoms whose basis functions are wanted without their nucleus. A
         !! group holding the attached end of a detached bond carries the
         !! detached atom this way: the functions describe the bond region, and
         !! the nucleus and its electrons stay with the fragment that owns them.

      if (present(ghost)) then
         call build_libcint_molecule(frag_z, frag_sym, frag_xyz, trim(opts%basis), mol, &
                                     error, ghost=ghost)
      else
         call build_libcint_molecule(frag_z, frag_sym, frag_xyz, trim(opts%basis), mol, error)
      end if
      if (error%has_error()) return
      call schwarz_bounds(mol, bounds, error)
   end subroutine open_fragment

   subroutine build_fragments(z, symbols, coords, owner, opts, frag, n_frag, afo, error, comm)
      !! Each fragment, and which of the others it needs the exact term for
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: owner(:)
      type(fmo_options_t), intent(in) :: opts
      type(fragment_t), allocatable, intent(out) :: frag(:)
      integer, intent(out) :: n_frag
      type(afo_context_t), intent(out) :: afo
      type(error_t), intent(inout) :: error
      type(comm_t), intent(in), optional :: comm

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

      if (opts%bond_breaking == "none") then
         call refuse_severed_bonds(z, coords, owner, n_atoms, error)
         if (error%has_error()) return
      else if (opts%bond_breaking == "afo") then
         ! Adjusted frozen orbitals. Restricted to the unembedded expansion for
         ! now, and that ordering is deliberate rather than incidental: a frozen
         ! orbital and an embedding field both describe the bond region, and the
         ! detached atom's share has to come out of the field before the two can
         ! be used together. Getting the fragmentation right first, where there
         ! is no field to double count against, is what the caps attempt should
         ! have done.
         if (opts%esp /= "none") then
            call error%set(ERROR_VALIDATION, "fmo: bond_breaking='afo' is implemented "// &
                           "for esp='none' only. A frozen orbital and an embedding "// &
                           "field both describe the detached bond, and removing the "// &
                           "detached atom's share of the field is not built yet -- it "// &
                           "is clean for point charges and not defined for an exact "// &
                           "density. Set keywords.fragmentation.embedding to 'none'")
            return
         end if
         call build_afo_context(z, symbols, coords, owner, opts, afo, error, comm)
         if (error%has_error()) return
      else
         call error%set(ERROR_VALIDATION, "fmo: bond_breaking='"// &
                        trim(opts%bond_breaking)//"' is not implemented for this "// &
                        "expansion. A cap has to be decided per n-mer rather than per "// &
                        "fragment -- a dimer spanning a cut bond must not carry the "// &
                        "caps its monomers needed. Use 'afo', which detaches a bond "// &
                        "with a frozen orbital instead, or the many-body path, which "// &
                        "caps covalent cuts already")
         return
      end if

      allocate (frag(n_frag))
      do f = 1, n_frag
         frag(f)%atoms = pack([(i, i=1, n_atoms)], owner == f)
         frag(f)%z = z(frag(f)%atoms)
         frag(f)%sym = symbols(frag(f)%atoms)
         frag(f)%xyz = coords(:, frag(f)%atoms)
         frag(f)%nelec = sum(frag(f)%z)
         ! A detached bond moves an electron between the two fragments it joins:
         ! the end that gets nothing of it is one short, the end that gets all
         ! of it is one over. Applied here so the closed-shell check below sees
         ! the count the fragment will actually be solved with -- ethane split
         ! into two methyls is 9 and 9 before this and 8 and 10 after it.
         if (afo%active) then
            frag(f)%nelec = frag(f)%nelec + group_electron_shift(afo%cuts, afo%n_cuts, [f])
         end if
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

   subroutine build_afo_context(z, symbols, coords, owner, opts, afo, error, comm)
      !! Every cut bond's frozen orbital, worked out once for the system
      !!
      !! A hybrid belongs to a bond and its surroundings, not to whoever is being
      !! solved, so this runs before any fragment does and costs one small SCF
      !! per cut bond however many fragments and n-mers there turn out to be.
      !!
      !! **Solved on one rank and shared, rather than everywhere.** The cost is
      !! not the reason -- a model system is a dozen atoms and its SCF is
      !! milliseconds, so recomputing it per rank wastes nothing anybody misses.
      !! The reason is that every rank has to end up with *bit-identical*
      !! hybrids. Each rank freezes orbitals in the fragments it owns and the
      !! energies are summed by an allreduce, so ranks that froze subtly
      !! different orbitals contribute from subtly different methods and the sum
      !! is quietly wrong -- no crash, nothing in the log.
      !!
      !! Recomputing would very probably agree, and "very probably" is the
      !! problem. Two things could separate ranks: unequal `OMP_NUM_THREADS`
      !! reassociates the BLAS reductions under the SCF, and the localization is
      !! a *discrete* choice on top of continuous data -- Jacobi sweeps take the
      !! largest pair gain, and near-degenerate gains on the symmetric little
      !! molecules model systems tend to be can flip on a rounding difference
      !! and give a different localized set rather than a perturbed one.
      !!
      !! **The failure path has to go through the collectives, not around
      !! them.** A leader that hits an error and returns leaves every other rank
      !! blocked forever at a reduction nobody reaches. So the leader records
      !! what went wrong as integers, every rank reduces them, and every rank
      !! reconstructs the same message from the result.
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: owner(:)
      type(fmo_options_t), intent(in) :: opts
      type(afo_context_t), intent(out) :: afo
      type(error_t), intent(inout) :: error
      type(comm_t), intent(in), optional :: comm

      type(system_geometry_t) :: geom
      type(afo_model_t) :: model
      type(afo_options_t) :: afo_opts
      type(error_t) :: local
      real(dp), allocatable :: hyb(:), flat(:)
      integer, allocatable :: lengths(:)
      integer :: status(3)
      integer :: i, n_on_bond, at, total

      afo%active = .false.
      if (error%has_error()) return

      geom%total_atoms = size(z)
      allocate (geom%element_numbers(size(z)), source=z)
      allocate (geom%coordinates(3, size(z)), source=coords)
      call find_severed_bonds(geom, owner, afo%cuts, afo%n_cuts)
      if (afo%n_cuts == 0) return

      do i = 1, afo%n_cuts
         if (afo%cuts(i)%in_ring) then
            call error%set(ERROR_VALIDATION, "fmo: the partition cuts a ring -- atoms "// &
                           to_char(afo%cuts(i)%atom_a)//" and "// &
                           to_char(afo%cuts(i)%atom_b)//" are still joined by another "// &
                           "path -- so the two fragments meet in more than one place. "// &
                           "A frozen orbital stands in for one detached bond and not "// &
                           "for a ring; fragment so that each pair of fragments is "// &
                           "joined at most once")
            return
         end if
         if (z(afo%cuts(i)%atom_a) == 1 .or. z(afo%cuts(i)%atom_b) == 1) then
            call error%set(ERROR_VALIDATION, "fmo: the partition detaches a bond at a "// &
                           "hydrogen (atoms "//to_char(afo%cuts(i)%atom_a)//" and "// &
                           to_char(afo%cuts(i)%atom_b)//"). A hydrogen has nothing left "// &
                           "to hybridise once its one bond is taken away; cut at a "// &
                           "heavy atom")
            return
         end if
      end do

      afo_opts%basis = opts%basis
      ! The rest of it, which was never copied. `afo_options_t` has carried
      ! these three fields all along and this line set only the basis, so the
      ! model system a cut bond's frozen orbital is localized from converged on
      ! that type's defaults no matter what the run was told.
      afo_opts%scf_max_iter = opts%scf_max_iter
      afo_opts%scf_energy_tol = opts%scf_energy_tol
      afo_opts%scf_density_tol = opts%scf_density_tol
      afo_opts%scf = opts%scf
      allocate (afo%hybrid(afo%n_cuts))
      allocate (afo%sym(size(symbols)), source=symbols)
      allocate (lengths(afo%n_cuts), source=0)
      status = 0

      ! `status` is [what went wrong, which bond, how many orbitals were on it],
      ! filled only by the rank that does the work and reduced below, so the
      ! message every rank raises is the same one.
      if (is_leader(comm)) then
         do i = 1, afo%n_cuts
            call build_afo_model(z, coords, afo%cuts(i), model, local)
            if (local%has_error()) then
               status = [1, i, 0]
               exit
            end if
            call bond_hybrid(model, afo_opts, hyb, n_on_bond, local)
            if (local%has_error()) then
               status = [2, i, 0]
               exit
            end if
            if (n_on_bond /= 1) then
               status = [3, i, n_on_bond]
               exit
            end if
            afo%hybrid(i)%coeff = hyb
            lengths(i) = size(hyb)
         end do
      end if

      if (spread_over(comm)) call allreduce(comm, status, 3, MPI_SUM)
      if (status(1) /= 0) then
         call afo_failure(afo%cuts(status(2)), status(1), status(3), local, error)
         return
      end if

      if (spread_over(comm)) then
         call allreduce(comm, lengths, afo%n_cuts, MPI_SUM)
         total = sum(lengths)
         allocate (flat(total), source=0.0_dp)
         if (is_leader(comm)) then
            at = 0
            do i = 1, afo%n_cuts
               flat(at + 1:at + lengths(i)) = afo%hybrid(i)%coeff
               at = at + lengths(i)
            end do
         end if
         call allreduce(comm, flat, total, MPI_SUM)
         at = 0
         do i = 1, afo%n_cuts
            if (allocated(afo%hybrid(i)%coeff)) deallocate (afo%hybrid(i)%coeff)
            allocate (afo%hybrid(i)%coeff(lengths(i)))
            afo%hybrid(i)%coeff = flat(at + 1:at + lengths(i))
            at = at + lengths(i)
         end do
         deallocate (flat)
      end if

      afo%active = .true.
      call logger%verbose("  fmo: "//to_char(afo%n_cuts)//" detached bond(s), one "// &
                          "frozen orbital each")
   end subroutine build_afo_context

   subroutine afo_failure(cut, kind, n_on_bond, local, error)
      !! The same message on every rank, rebuilt from the reduced status
      !!
      !! `local` carries the lead rank's own error text and is the better
      !! message where it exists, but it exists only on that rank. So the text
      !! is reconstructed from integers that every rank has, and the leader's
      !! own message is appended where it has one -- identical failures
      !! everywhere, with detail where detail is available.
      type(severed_bond_t), intent(in) :: cut
      integer, intent(in) :: kind, n_on_bond
      type(error_t), intent(inout) :: local
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: bond

      bond = "the bond between atoms "//to_char(cut%atom_a)//" and "// &
             to_char(cut%atom_b)

      select case (kind)
      case (1)
         call error%set(ERROR_VALIDATION, "fmo: the model system for "//bond// &
                        " could not be built")
      case (2)
         call error%set(ERROR_VALIDATION, "fmo: the model system for "//bond// &
                        " could not be solved, so there is no orbital to freeze")
      case default
         call error%set(ERROR_VALIDATION, "fmo: "//to_char(n_on_bond)//" localized "// &
                        "orbitals sit on "//bond//", so it is not a single bond. One "// &
                        "frozen orbital stands in for one electron pair; cut at a "// &
                        "single bond")
      end select
      if (local%has_error()) then
         call logger%verbose("  fmo: "//trim(local%get_message()))
      end if
   end subroutine afo_failure

   subroutine assemble_group(frag, members, afo, z, coords, group, error)
      !! One n-mer's geometry, electron count and boundaries
      !!
      !! The members end to end, then any ghosts they need. Ghosts go last on
      !! purpose: every member keeps the contiguous run of basis functions the
      !! block bookkeeping everywhere else assumes, and what is new sits past
      !! the end of it.
      !!
      !! The boundary set comes from `cuts_outside_group`, computed from this
      !! group's own members every time. A bond cut between two monomers is
      !! whole inside the dimer that holds both, so that dimer has no boundary
      !! there -- and inheriting one from its members is the mistake that cost
      !! 11 Hartree when this was tried with caps.
      type(fragment_t), intent(in) :: frag(:)
      integer, intent(in) :: members(:)
      type(afo_context_t), intent(in) :: afo
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords(:, :)
      type(group_t), intent(out) :: group
      type(error_t), intent(inout) :: error

      integer, allocatable :: outside(:), slot_of(:)
      logical, allocatable :: ghosted(:)
      integer :: n_outside, n_ghost, n_atoms, m, at, n_here, i, c, bda, slot

      if (error%has_error()) return

      group%n_real = 0
      do m = 1, size(members)
         group%n_real = group%n_real + size(frag(members(m))%z)
      end do

      n_outside = 0
      n_ghost = 0
      allocate (ghosted(size(z)), source=.false.)
      if (afo%active) then
         call cuts_outside_group(afo%cuts, afo%n_cuts, members, outside, n_outside)
         do i = 1, n_outside
            c = outside(i)
            ! Holding the attached end means the detached atom is somebody
            ! else's, so its functions have to be brought in without its nucleus.
            if (any(members == afo%cuts(c)%frag_a)) cycle
            ! Counted once per *atom*, not once per boundary. One atom can be
            ! the detached end of two bonds -- a middle carbon numbered below
            ! both its neighbours -- and a group holding the far end of both
            ! would otherwise bring it in twice. Two copies of one atom's
            ! functions make the overlap exactly singular; the canonical
            ! orthogonalisation absorbs that, so the answer survives, but it
            ! survives by leaning on a linear-dependence threshold instead of
            ! by being right.
            if (ghosted(afo%cuts(c)%atom_a)) cycle
            ghosted(afo%cuts(c)%atom_a) = .true.
            n_ghost = n_ghost + 1
         end do
         ghosted = .false.
      end if

      n_atoms = group%n_real + n_ghost
      allocate (group%z(n_atoms), group%sym(n_atoms), group%xyz(3, n_atoms))
      allocate (group%ghost(n_atoms), source=.false.)
      allocate (slot_of(size(z)), source=0)

      at = 0
      do m = 1, size(members)
         n_here = size(frag(members(m))%z)
         group%z(at + 1:at + n_here) = frag(members(m))%z
         group%sym(at + 1:at + n_here) = frag(members(m))%sym
         group%xyz(:, at + 1:at + n_here) = frag(members(m))%xyz
         do i = 1, n_here
            slot_of(frag(members(m))%atoms(i)) = at + i
         end do
         at = at + n_here
      end do

      group%nelec = sum(group%z(:group%n_real))
      group%n_bound = n_outside
      allocate (group%bda_slot(max(n_outside, 1)))
      allocate (group%occupied(max(n_outside, 1)))
      allocate (group%cut_of(max(n_outside, 1)))

      do i = 1, n_outside
         c = outside(i)
         bda = afo%cuts(c)%atom_a
         group%cut_of(i) = c
         if (any(members == afo%cuts(c)%frag_a)) then
            ! We own the detached atom and get nothing of the bond: its hybrid
            ! is frozen empty.
            group%occupied(i) = .false.
            slot = slot_of(bda)
            if (slot == 0) then
               call error%set(ERROR_VALIDATION, "fmo: a detached atom is not among the "// &
                              "atoms of the fragment that owns it")
               return
            end if
         else
            ! We own the attached end and get all of the bond: bring the
            ! detached atom's functions in as a ghost and hold its hybrid full.
            ! One copy per atom -- a second boundary on the same detached atom
            ! puts its hybrid on the copy already there, where the two are
            ! different orbitals of one atom and independent for that reason.
            group%occupied(i) = .true.
            if (ghosted(bda)) then
               slot = slot_of(bda)
            else
               at = at + 1
               slot = at
               group%z(slot) = z(bda)
               group%sym(slot) = afo%sym(bda)
               group%xyz(:, slot) = coords(:, bda)
               group%ghost(slot) = .true.
               ghosted(bda) = .true.
               slot_of(bda) = slot
            end if
         end if
         group%bda_slot(i) = slot
      end do

      if (afo%active) then
         group%nelec = group%nelec + group_electron_shift(afo%cuts, afo%n_cuts, members)
      end if

      deallocate (slot_of, ghosted)
   end subroutine assemble_group

   subroutine group_projector(group, mol, afo, proj, active, error)
      !! The constraint this group's boundaries imply, if it has any
      type(group_t), intent(in) :: group
      type(libcint_molecule_t), intent(in) :: mol
      type(afo_context_t), intent(in) :: afo
      type(fock_projector_t), intent(out) :: proj
      logical, intent(out) :: active
      type(error_t), intent(inout) :: error

      type(afo_hybrid_t), allocatable :: hybrids(:)
      real(dp), allocatable :: frozen(:, :), basis(:, :), s(:, :)
      integer :: i, n_frozen_occ, n_mo

      active = .false.
      if (error%has_error()) return
      if (group%n_bound == 0) return

      allocate (hybrids(group%n_bound))
      do i = 1, group%n_bound
         hybrids(i)%coeff = afo%hybrid(group%cut_of(i))%coeff
      end do

      call build_group_frozen(mol, group%bda_slot(:group%n_bound), &
                              group%occupied(:group%n_bound), hybrids, frozen, &
                              n_frozen_occ, error)
      if (error%has_error()) return

      call mol%overlap(s)
      call build_frozen_basis(frozen, n_frozen_occ, s, basis, n_mo, error)
      if (error%has_error()) return
      call proj%init(basis, s, n_frozen_occ, group%n_bound, AFO_SHIFT, error)
      if (error%has_error()) return
      active = .true.
   end subroutine group_projector

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
         ! Only the real atoms map back. A cap's charge belongs to no atom
         ! of the system, and `charges` is as long as the molecule the SCF
         ! saw, which includes them.
         q_all(frag(f)%atoms) = frag(f)%charges(1:size(frag(f)%atoms))
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

   subroutine nmer_term(frag, n_frag, members, z, coords, q_all, opts, afo, &
                        e_internal, e_resp, all_converged, error)
      !! One n-mer, in the field of every fragment outside it
      !!
      !! The pair case generalised: nothing here knew that a group was two
      !! fragments except the block bookkeeping, and that was always a running
      !! offset rather than a front half and a back half.
      type(fragment_t), intent(in) :: frag(:)
      integer, intent(in) :: n_frag
      integer, intent(in) :: members(:)
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords(:, :)
      real(dp), allocatable, intent(in) :: q_all(:)
      type(fmo_options_t), intent(in) :: opts
      type(afo_context_t), intent(in) :: afo
      real(dp), intent(out) :: e_internal, e_resp
      logical, intent(inout) :: all_converged
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      integer, allocatable :: deck_guess
      type(group_t) :: group
      type(fock_projector_t) :: proj
      real(dp), allocatable :: bounds(:, :), d_split(:, :), u(:, :)
      logical, allocatable :: inside(:)
      integer, allocatable :: near(:), ao_off(:), ao_count(:)
      integer :: m, at, nao_m, expect, nelec
      logical :: held

      ! The n-mer's geometry, its fragments end to end in the order given, and
      ! then whatever ghosts its own boundaries call for.
      expect = 0
      do m = 1, size(members)
         expect = expect + frag(members(m))%nao
      end do

      call assemble_group(frag, members, afo, z, coords, group, error)
      if (error%has_error()) return
      nelec = group%nelec

      call open_fragment(group%z, group%sym, group%xyz, opts, mol, bounds, error, &
                         ghost=group%ghost)
      if (error%has_error()) return

      ! libcint orders functions by atom and the atoms went in fragment by
      ! fragment, so each member owns a contiguous run. Worth checking rather
      ! than assuming: a silent mismatch would be a plausible wrong answer
      ! rather than a failure.
      ! The members' functions must still be the leading block, whatever ghosts
      ! follow them. Counted over the real atoms rather than compared against
      ! `mol%nao`, which now includes the ghosts.
      allocate (ao_off(mol%natm), ao_count(mol%natm))
      call atom_ao_blocks(mol, ao_off, ao_count)
      if (sum(ao_count(:group%n_real)) /= expect) then
         call error%set(ERROR_VALIDATION, "fmo: the n-mer basis is not its fragments' "// &
                        "bases end to end ("//to_char(expect)//" /= "// &
                        to_char(sum(ao_count(:group%n_real)))//")")
         return
      end if
      deallocate (ao_off, ao_count)

      ! The member densities side by side: what the n-mer's own Coulomb
      ! contribution is removed with, and what dD is measured against.
      allocate (d_split(mol%nao, mol%nao), source=0.0_dp)
      at = 0
      do m = 1, size(members)
         nao_m = frag(members(m))%nao
         d_split(at + 1:at + nao_m, at + 1:at + nao_m) = frag(members(m))%density
         at = at + nao_m
      end do

      allocate (inside(size(z)), source=.false.)
      do m = 1, size(members)
         inside(frag(members(m))%atoms) = .true.
      end do

      call near_fragments(frag, n_frag, members, effective_resppc(opts), near, error)
      if (error%has_error()) return

      call embedding_operator(mol, group%z, group%sym, group%xyz, near, frag, n_frag, &
                              inside, z, coords, q_all, opts, u, error)
      if (error%has_error()) return

      call group_projector(group, mol, afo, proj, held, error)
      if (error%has_error()) return

      ! Resolved before the branch, not inside it. `deck_guess` is passed on all
      ! four, but only one used to fill it -- and an unallocated allocatable
      ! arrives at an optional dummy as *absent*, so `keywords.scf.guess` was
      ! silently dropped on every path but the embedded-and-constrained one.
      ! `fmo_guess_kind` leaves it unallocated when the deck said nothing or
      ! said `auto`, which is the intended "absent", so hoisting changes only
      ! the decks that asked for a specific guess.
      call fmo_guess_kind(opts, deck_guess)

      if (allocated(u) .and. held) then
         call run_libcint_rhf(mol, nelec, opts%scf_max_iter, opts%scf_energy_tol, &
                              opts%scf_density_tol, show_inner_scf(), scf, error, scf=opts%scf, guess=deck_guess, &
                              h_extra=u, projector=proj)
      else if (allocated(u)) then
         call run_libcint_rhf(mol, nelec, opts%scf_max_iter, opts%scf_energy_tol, &
                            opts%scf_density_tol, show_inner_scf(), scf, error, scf=opts%scf, guess=deck_guess, h_extra=u)
      else if (held) then
         call run_libcint_rhf(mol, nelec, opts%scf_max_iter, opts%scf_energy_tol, &
                              opts%scf_density_tol, show_inner_scf(), scf, error, scf=opts%scf, guess=deck_guess, &
                              projector=proj)
      else
         call run_libcint_rhf(mol, nelec, opts%scf_max_iter, opts%scf_energy_tol, &
                              opts%scf_density_tol, show_inner_scf(), scf, error, scf=opts%scf, guess=deck_guess)
      end if
      if (error%has_error()) return
      if (.not. scf%converged) all_converged = .false.

      e_internal = scf%energy
      e_resp = 0.0_dp
      if (allocated(u)) then
         if (opts%expansion /= "mbe") then
            e_internal = e_internal - sum(scf%density*u)
            e_resp = sum((scf%density - d_split)*u)
         end if
      end if
   end subroutine nmer_term

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

   subroutine solve_fragment(frag, n_frag, which, z, coords, q_all, opts, afo, &
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
      type(afo_context_t), intent(in) :: afo
      logical, intent(inout) :: all_converged
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: bare
         !! Solve it in vacuum. True only for the very first pass, which has no
         !! field to be solved against because no fragment has a density yet --
         !! the field is what the passes after it are for.

      type(libcint_molecule_t) :: mol
      type(group_t) :: group
      type(fock_projector_t) :: proj
      real(dp), allocatable :: bounds(:, :), u(:, :), q(:)
      logical, allocatable :: inside(:)
      logical :: isolated, held

      isolated = .false.
      if (present(bare)) isolated = bare

      ! A monomer is a group of one, and goes through the same assembly as an
      ! n-mer so that its boundaries are decided the same way.
      call assemble_group(frag, [which], afo, z, coords, group, error)
      if (error%has_error()) return

      call open_fragment(group%z, group%sym, group%xyz, opts, mol, bounds, error, &
                         ghost=group%ghost)
      if (error%has_error()) return

      call group_projector(group, mol, afo, proj, held, error)
      if (error%has_error()) return

      if (.not. isolated) then
         allocate (inside(size(z)), source=.false.)
         inside(frag(which)%atoms) = .true.
         call embedding_operator(mol, frag(which)%z, frag(which)%sym, frag(which)%xyz, &
                                 frag(which)%near, frag, n_frag, inside, z, coords, &
                                 q_all, opts, u, error)
         if (error%has_error()) return
      end if

      call inner_scf(frag(which), mol, opts, error, u, all_converged, proj, held)
      if (error%has_error()) return

      ! A fragment carrying ghosts was solved in a bigger basis than it owns,
      ! and everything downstream indexes its density by `nao`, which counts
      ! only its own atoms. Ghosts are appended last, so the block it owns is
      ! the leading corner; keeping that keeps `d_split` and the block layout
      ! valid. The ghost block is dropped rather than stored, which is sound
      ! while this is restricted to `esp = "none"` -- nothing reads a monomer
      ! density there. An embedded version will want the whole thing, and will
      ! have to widen the layout rather than truncate here.
      if (size(frag(which)%density, 1) /= frag(which)%nao) then
         frag(which)%density = frag(which)%density(:frag(which)%nao, :frag(which)%nao)
      end if

      if (opts%esp /= "none") then
         call fragment_charges(mol, frag(which)%density, opts%far_field, q, error)
         if (error%has_error()) return
         frag(which)%charges = q
      end if
   end subroutine solve_fragment

   subroutine calculate_monomers(frag, n_frag, z, coords, opts, afo, res, all_converged, error, comm)
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
      type(afo_context_t), intent(in) :: afo
      type(fmo_result_t), intent(inout) :: res
      logical, intent(inout) :: all_converged
      type(error_t), intent(inout) :: error
      type(comm_t), intent(in), optional :: comm

      real(dp), allocatable :: q_all(:)
      real(dp) :: e_sum, e_prev
      integer :: i, outer
      logical :: show_conv

      ! Isolated fragments, to start the field from something.
      call all_charges(frag, n_frag, size(z), opts, q_all, error)
      if (error%has_error()) return
      do i = 1, n_frag
         if (.not. mine(i, comm)) cycle
         call solve_fragment(frag, n_frag, i, z, coords, q_all, opts, afo, all_converged, &
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
      show_conv = show_outer(comm)
      call convergence_header(show_conv, "FMO monomer SCF", &
                              "    iter            monomer sum        change", 45)
      do outer = 1, opts%max_outer
         call all_charges(frag, n_frag, size(z), opts, q_all, error)
         if (error%has_error()) return

         ! Independent within a pass: every fragment reads the densities the
         ! last pass left and none reads this pass's, so who computes which is
         ! free. The exchange after is the barrier.
         do i = 1, n_frag
            if (.not. mine(i, comm)) cycle
            call solve_fragment(frag, n_frag, i, z, coords, q_all, opts, afo, &
                                all_converged, error)
            if (error%has_error()) return
         end do
         call exchange_monomers(frag, n_frag, comm)

         e_sum = sum(frag(:)%energy)
         res%outer_iterations = outer
         res%outer_change = abs(e_sum - e_prev)
         call fmo_outer_row(show_conv, outer, e_sum, res%outer_change)
         if (res%outer_change < opts%outer_tol) then
            res%converged = .true.
            call convergence_footer(show_conv, .true., outer, "outer iterations", 45)
            return
         end if
         e_prev = e_sum
      end do

      call convergence_footer(show_conv, .false., opts%max_outer, "outer iterations", 45)
      call error%set(ERROR_VALIDATION, "fmo: the outer SCF did not settle in "// &
                     to_char(opts%max_outer)//" passes; the monomer sum was still "// &
                     "moving by "//to_char(res%outer_change)//" Hartree")
   end subroutine calculate_monomers

   subroutine calculate_polymers(frag, n_frag, z, coords, opts, afo, res, all_converged, error, comm)
      !! Every n-mer from pairs up to the truncation level
      !!
      !! Independent of each other and of everything else once the monomers have
      !! settled, so one bag of tasks with no barrier inside it -- and the
      !! expensive part, since the count is C(N,n) rather than N.
      !!
      !! **The correction each n-mer contributes** is its internal energy less
      !! everything its own subsets already accounted for:
      !!
      !!     dE_S = E'_S - sum over proper non-empty subsets T of S of dE_T
      !!
      !! which unrolls to `E'_IJ - E'_I - E'_J` for a pair and to the usual
      !! three-body expression for a trimer, without either being written down.
      !! Summing dE over every subset up to the level *is* the truncated
      !! expansion, so the total is one sum and not a per-level formula. When
      !! the level reaches the fragment count the corrections telescope and the
      !! result is the supermolecular energy exactly.
      !!
      !! The response term of an n-mer is part of that n-mer's correction and so
      !! sits inside the recursion. Left outside it, a pair's response would
      !! never be cancelled by the trimer containing it, and the total would
      !! miss exactness by exactly the sum of them -- which is how this was
      !! found. `response_sum` is reported separately but is not added again.
      type(fragment_t), intent(inout) :: frag(:)
      integer, intent(in) :: n_frag
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords(:, :)
      type(fmo_options_t), intent(in) :: opts
      type(afo_context_t), intent(in) :: afo
      type(fmo_result_t), intent(inout) :: res
      logical, intent(inout) :: all_converged
      type(error_t), intent(inout) :: error
      type(comm_t), intent(in), optional :: comm

      real(dp), allocatable :: q_all(:), totals(:)
      integer, allocatable :: terms(:, :), term_size(:)
      real(dp), allocatable :: correction(:)
      real(dp) :: e_internal, e_resp
      integer :: n_terms, t, task, level, n_nmers

      level = min(opts%level, n_frag)
      if (level < 2) return

      call all_charges(frag, n_frag, size(z), opts, q_all, error)
      if (error%has_error()) return

      call enumerate_terms(n_frag, level, terms, term_size, n_terms)
      allocate (correction(n_terms), source=0.0_dp)

      ! Count the n-mers (size >= 2) so the progress below has a denominator: the
      ! monomers are already solved, so they are not part of this phase's work.
      n_nmers = count(term_size(1:n_terms) >= 2)
      if (is_leader(comm)) then
         call logger%info("  fmo: "//to_char(n_nmers)//" n-mers up to level "// &
                          to_char(level))
      end if

      ! Monomers are level one and already solved; their correction is just
      ! their energy, which is what the recursion below subtracts against.
      do t = 1, n_terms
         if (term_size(t) /= 1) cycle
         if (opts%expansion == "mbe") then
            correction(t) = frag(terms(1, t))%energy_total
         else
            correction(t) = frag(terms(1, t))%energy
         end if
      end do

      task = 0
      do t = 1, n_terms
         if (term_size(t) < 2) cycle
         task = task + 1

         ! Coarse progress at info level so a long n-mer sweep shows movement
         ! without the per-SCF detail. On the loop index, not the owned subset, so
         ! it advances smoothly; leader-guarded so it prints once under MPI.
         if (is_leader(comm)) then
            if (mod(task, max(1, n_nmers/10)) == 0 .or. task == n_nmers) then
               call logger%info("  fmo: n-mer "//to_char(task)//"/"//to_char(n_nmers))
            end if
         end if

         if (.not. mine(task, comm)) cycle

         call nmer_term(frag, n_frag, terms(1:term_size(t), t), z, coords, q_all, &
                        opts, afo, e_internal, e_resp, all_converged, error)
         if (error%has_error()) return
         ! The response goes inside the correction, not alongside it. A larger
         ! n-mer subtracts its subsets' corrections whole, so a response left
         ! outside the recursion never gets cancelled and survives into a total
         ! that should have telescoped exactly. At level two nothing contains a
         ! pair so the two placements agree, which is why this only shows up
         ! once there are trimers.
         correction(t) = e_internal + e_resp
         res%response_sum = res%response_sum + e_resp
      end do

      if (spread_over(comm)) then
         ! Each rank filled only its own n-mers, the rest left at zero, so a sum
         ! gathers them. Monomer entries were filled by every rank, so they are
         ! zeroed on all but one first or they would be counted size() times.
         if (comm%rank() /= 0) then
            do t = 1, n_terms
               if (term_size(t) == 1) correction(t) = 0.0_dp
            end do
         end if
         allocate (totals(n_terms + 1))
         totals(1:n_terms) = correction
         totals(n_terms + 1) = res%response_sum
         call allreduce(comm, totals, size(totals), MPI_SUM)
         correction = totals(1:n_terms)
         res%response_sum = totals(n_terms + 1)
      end if

      ! Subtract what the subsets already covered. Ordered by size, so every
      ! subset of a term is final before the term is reduced.
      call subtract_subsets(terms, term_size, n_terms, correction)

      res%pair_sum = 0.0_dp
      do t = 1, n_terms
         if (term_size(t) >= 2) res%pair_sum = res%pair_sum + correction(t)
      end do
   end subroutine calculate_polymers

   subroutine enumerate_terms(n_frag, level, terms, term_size, n_terms)
      !! Every combination of fragments from one up to `level`, smallest first
      !!
      !! Size order matters: the correction for a term is reduced by its
      !! subsets, so each subset has to be final before anything containing it
      !! is touched.
      integer, intent(in) :: n_frag, level
      integer, allocatable, intent(out) :: terms(:, :), term_size(:)
      integer, intent(out) :: n_terms

      integer, allocatable :: pick(:)
      integer :: m, total, k

      total = 0
      do m = 1, level
         total = total + n_choose(n_frag, m)
      end do
      allocate (terms(level, total), source=0)
      allocate (term_size(total), source=0)

      n_terms = 0
      do m = 1, level
         allocate (pick(m))
         do k = 1, m
            pick(k) = k
         end do
         do
            n_terms = n_terms + 1
            terms(1:m, n_terms) = pick
            term_size(n_terms) = m
            if (.not. step_combination(pick, m, n_frag)) exit
         end do
         deallocate (pick)
      end do
   end subroutine enumerate_terms

   function step_combination(pick, m, n) result(more)
      !! Advance a combination in lexicographic order; false when exhausted
      integer, intent(inout) :: pick(:)
      integer, intent(in) :: m, n
      logical :: more

      integer :: i, k

      more = .false.
      do i = m, 1, -1
         if (pick(i) < n - m + i) then
            pick(i) = pick(i) + 1
            do k = i + 1, m
               pick(k) = pick(k - 1) + 1
            end do
            more = .true.
            return
         end if
      end do
   end function step_combination

   function n_choose(n, k) result(c)
      !! Binomial coefficient, built up rather than from factorials
      integer, intent(in) :: n, k
      integer :: c

      integer :: i

      c = 1
      do i = 1, k
         c = c*(n - k + i)/i
      end do
   end function n_choose

   subroutine subtract_subsets(terms, term_size, n_terms, correction)
      !! Reduce each term by the corrections its own subsets already carry
      integer, intent(in) :: terms(:, :)
      integer, intent(in) :: term_size(:)
      integer, intent(in) :: n_terms
      real(dp), intent(inout) :: correction(:)

      integer :: t, u

      do t = 1, n_terms
         if (term_size(t) < 2) cycle
         do u = 1, n_terms
            if (term_size(u) >= term_size(t)) cycle
            if (.not. is_subset(terms(1:term_size(u), u), terms(1:term_size(t), t))) cycle
            correction(t) = correction(t) - correction(u)
         end do
      end do
   end subroutine subtract_subsets

   function is_subset(small, big) result(inside)
      !! Whether every member of `small` appears in `big`
      integer, intent(in) :: small(:), big(:)
      logical :: inside

      integer :: i

      inside = .true.
      do i = 1, size(small)
         if (.not. any(big == small(i))) then
            inside = .false.
            return
         end if
      end do
   end function is_subset

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

   function is_leader(comm) result(leads)
      !! Whether this rank should emit the run's shared log lines
      !!
      !! The outer-convergence and progress lines report reduced totals, so they
      !! are identical on every rank -- without this each rank prints its own copy
      !! and the log comes back interleaved N times under MPI. The per-fragment
      !! SCF tables are deliberately *not* gated this way: each rank solves a
      !! different set of fragments, so letting every rank print its own is the
      !! point. Absent `comm` is a single rank, which always leads.
      type(comm_t), intent(in), optional :: comm
      logical :: leads

      leads = .true.
      if (.not. spread_over(comm)) return
      leads = comm%rank() == 0
   end function is_leader

   function show_inner_scf() result(show)
      !! Whether to print each fragment's per-iteration SCF table
      !!
      !! The individual monomer and n-mer SCFs are the finest thing this method
      !! prints and there are a great many of them, so they sit at the deepest
      !! level. Gating on the logger level rather than a threaded flag keeps the
      !! decision in one place and matches how the MBE path decides the same.
      !TODO: replace debug_level with verbose_level once pic defines large_info;
      !      the outer convergence (see calculate_monomers) drops to large_info at
      !      the same time, so inner SCF and outer stay one level apart.
      logical :: show
      integer :: level

      call logger%configuration(level=level)
      show = level >= debug_level
   end function show_inner_scf

   function show_outer(comm) result(show)
      !! Whether to print the FMO outer-loop (monomer SCF) convergence table
      !!
      !! Leader-guarded -- the monomer sum is a reduced total, identical on every
      !! rank -- and gated at verbose: one level above the info progress lines and
      !! one below the debug per-fragment SCF tables. The level travels in this
      !! flag rather than in the logger call because the shared table frame
      !! (convergence_header/footer) prints through logger%info for every table.
      !TODO: replace verbose_level with large_info once pic defines it, so the
      !      outer table sits just below the inner SCF tables (see show_inner_scf).
      type(comm_t), intent(in), optional :: comm
      logical :: show
      integer :: level

      call logger%configuration(level=level)
      show = is_leader(comm) .and. level >= verbose_level
   end function show_outer

   subroutine fmo_outer_row(show, iter, monomer_sum, change)
      !! One outer iteration's line: the monomer energy sum and how far it moved
      !!
      !! The FMO loop's row for the shared convergence table -- the SCF row minus
      !! the DIIS depth and per-iteration timings an outer loop has no analogue
      !! for. The frame around it is convergence_header/convergence_footer.
      logical, intent(in) :: show
      integer, intent(in) :: iter
      real(dp), intent(in) :: monomer_sum, change

      character(len=100) :: line

      if (.not. show) return
      write (line, "(i8,f23.12,es14.3)") iter, monomer_sum, change
      call logger%info(trim(line))
   end subroutine fmo_outer_row

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
               ! Sliced, as `all_charges` slices. A fragment solved with ghosts
               ! has a charge per atom of the molecule it saw, which is more
               ! than the atoms it owns, and the buffer is laid out by the
               ! latter. Unreachable while a detached bond forces `esp="none"`
               ! and no charges are computed at all, and left defensive because
               ! the assignment would otherwise be a shape mismatch the moment
               ! that restriction lifts.
               buf(at + n + 3:at + n + 2 + size(frag(f)%atoms)) = &
                  frag(f)%charges(1:size(frag(f)%atoms))
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
         frag(f)%charges(1:size(frag(f)%atoms)) = &
            buf(at + n + 3:at + n + 2 + size(frag(f)%atoms))
         at = at + n + 2 + size(frag(f)%atoms)
      end do
   end subroutine exchange_monomers

   subroutine inner_scf(f, mol, opts, error, u, all_converged, proj, held)
      !! The inner SCF: this fragment's orbitals, against a fixed external field
      type(fragment_t), intent(inout) :: f
      type(libcint_molecule_t), intent(in) :: mol
      type(fmo_options_t), intent(in) :: opts
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(in), optional :: u(:, :)
      logical, intent(inout), optional :: all_converged
      type(fock_projector_t), intent(in), optional :: proj
         !! The boundary constraint, when this fragment has a detached bond
      logical, intent(in), optional :: held
         !! Whether `proj` carries anything. Separate from its presence because
         !! a fragment with no boundary still has a projector object.

      type(rhf_result_t) :: scf
      integer, allocatable :: deck_guess
      logical :: embedded, constrained

      embedded = .false.
      if (present(u)) embedded = allocated(u)
      constrained = .false.
      if (present(held)) constrained = held

      ! Resolved before the branch, not inside it. `deck_guess` is passed on all
      ! four, but only one used to fill it -- and an unallocated allocatable
      ! arrives at an optional dummy as *absent*, so `keywords.scf.guess` was
      ! silently dropped on every path but the embedded-and-constrained one.
      ! `fmo_guess_kind` leaves it unallocated when the deck said nothing or
      ! said `auto`, which is the intended "absent", so hoisting changes only
      ! the decks that asked for a specific guess.
      call fmo_guess_kind(opts, deck_guess)

      if (embedded .and. constrained) then
         call run_libcint_rhf(mol, f%nelec, opts%scf_max_iter, opts%scf_energy_tol, &
                              opts%scf_density_tol, show_inner_scf(), scf, error, scf=opts%scf, guess=deck_guess, &
                              h_extra=u, projector=proj)
      else if (embedded) then
         call run_libcint_rhf(mol, f%nelec, opts%scf_max_iter, opts%scf_energy_tol, &
                            opts%scf_density_tol, show_inner_scf(), scf, error, scf=opts%scf, guess=deck_guess, h_extra=u)
      else if (constrained) then
         call run_libcint_rhf(mol, f%nelec, opts%scf_max_iter, opts%scf_energy_tol, &
                              opts%scf_density_tol, show_inner_scf(), scf, error, scf=opts%scf, guess=deck_guess, &
                              projector=proj)
      else
         call run_libcint_rhf(mol, f%nelec, opts%scf_max_iter, opts%scf_energy_tol, &
                              opts%scf_density_tol, show_inner_scf(), scf, error, scf=opts%scf, guess=deck_guess)
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
