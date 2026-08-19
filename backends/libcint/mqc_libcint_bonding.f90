!! Bonding analysis on a converged wave function
module mqc_libcint_bonding
   !! One entry point that takes a converged set of molecular orbitals and
   !! produces the quasi-atomic bonding picture: atomic populations and charges,
   !! then every orbital pair worth reporting with its bond order, its kinetic
   !! bond order and what kind of orbital each end is.
   !!
   !! The pipeline itself lives in `mqc_aambs`, `mqc_libcint_quao` and
   !! `mqc_libcint_quao_report`, each validated on its own. What is here is only
   !! the assembly -- which spaces to build, in what order, and what to hand to
   !! the report -- so that a deck can ask for the analysis without a caller
   !! having to know that the valence-internal space is the occupied valence
   !! orbitals followed by the valence-virtual ones.
   !!
   !! **This is a Hartree-Fock analysis for now.** The construction takes a
   !! one-particle density matrix and does not care where it came from, and
   !! `quasi_atomic_orbitals` accepts one explicitly for exactly that reason.
   !! But nothing here supplies a correlated density yet, so what is analysed is
   !! the reference determinant.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_io, only: to_char
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_aambs, only: aambs_dimensions, aambs_dimensions_t, aambs_element_counts
   use mqc_libcint_integrals, only: libcint_molecule_t, mixed_basis_overlap
   use mqc_libcint_quao, only: build_aambs_molecule, mo_aambs_overlap, &
                               valence_virtual_orbitals, vvo_result_t, &
                               aambs_atom_ranges, quasi_atomic_orbitals, &
                               quao_result_t, orient_quasi_atomic_orbitals, &
                               kinetic_bond_orders
   use mqc_physical_constants, only: HARTREE_TO_KCALMOL
   use mqc_libcint_atomic_guess, only: free_atom_energies
   use mqc_libcint_ieda, only: kinetic_decomposition, print_kinetic_decomposition, &
                               nuclear_attraction_per_atom, quao_nuclear_attraction, &
                               nuclear_decomposition, print_nuclear_decomposition, &
                               combine_quao_sets, quao_interference, kinetic_total, &
                               quao_eris, two_electron_decomposition, &
                               print_two_electron_decomposition, &
                               nuclear_repulsion_pairs, print_total_decomposition, &
                               print_interatomic_split, &
                               active_cumulant, quao_projection, transform_cumulant
   use mqc_libcint_quao_report, only: quao_labels_t, label_quasi_atomic_orbitals, &
                                      print_quao_report
   implicit none
   private

   public :: run_quao_analysis
   public :: BONDING_NONE, BONDING_GMS_QUAO
   public :: bonding_analysis_kind
   public :: bonding_analysis_name

   integer, parameter :: BONDING_NONE = 0
   integer, parameter :: BONDING_GMS_QUAO = 1
      !! The Ruedenberg quasi-atomic analysis as GAMESS implements it. Named for
      !! the reference implementation rather than the papers because the labels
      !! and the thresholds are GAMESS's, and only the underlying quantities are
      !! the papers'.

contains

   pure function bonding_analysis_kind(name) result(kind)
      !! Parse a deck's `properties.bonding_analysis` value
      character(len=*), intent(in) :: name
      integer :: kind

      select case (trim(adjustl(name)))
      case ("", "none")
         kind = BONDING_NONE
      case ("gms_quao", "quao")
         kind = BONDING_GMS_QUAO
      case default
         kind = -1
      end select
   end function bonding_analysis_kind

   pure function bonding_analysis_name(kind) result(name)
      !! The spelling a deck would use, for messages
      integer, intent(in) :: kind
      character(len=:), allocatable :: name

      select case (kind)
      case (BONDING_GMS_QUAO)
         name = "gms_quao"
      case default
         name = "none"
      end select
   end function bonding_analysis_name

   subroutine run_quao_analysis(mol, atomic_numbers, element_symbols, coordinates, &
                                orbitals, n_electrons, error, verbose, threshold, &
                                occupations, active_orbitals, active_dm1, active_dm2, &
                                reference_energy)
      !! The quasi-atomic bonding analysis, start to finish
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)     !! (3, natm), Bohr
      real(dp), intent(in) :: orbitals(:, :)        !! Converged, (n_ao, n_mo)
      integer, intent(in) :: n_electrons
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: verbose
      real(dp), intent(in), optional :: threshold
         !! kcal/mol. Pairs whose kinetic bond order is weaker are counted and
         !! not printed.
      real(dp), intent(in), optional :: occupations(:)
         !! (n_mo) occupation of each orbital, which must then be a *natural*
         !! orbital -- the density has to be diagonal in the basis given for
         !! these numbers to be its whole content.
         !!
         !! Absent means the reference determinant: two electrons in each of the
         !! first `n_occupied` orbitals. Present is what a correlated wave
         !! function supplies, and is the only difference between analysing one
         !! and analysing the other.

      type(libcint_molecule_t) :: aambs
      type(aambs_dimensions_t) :: dims
      type(vvo_result_t) :: vvo
      type(quao_result_t) :: quao
      type(quao_labels_t) :: labels
      real(dp), intent(in), optional :: active_orbitals(:, :)
         !! (n_ao, n_active), the orbitals `active_dm2` is expressed in. These
         !! are the optimised active orbitals, not the natural ones the rest of
         !! the analysis runs on: the two-particle density belongs to the basis
         !! it was computed in, and is projected from there.
      real(dp), intent(in), optional :: active_dm1(:, :)
      real(dp), intent(in), optional :: active_dm2(:, :, :, :)
         !! Present together with the two above, or not at all. Their presence
         !! is what makes the energy decomposition correlated rather than a
         !! single-determinant one.
      real(dp), intent(in), optional :: reference_energy
         !! The energy the calculation reported, to check the decomposition
         !! against. Required for a correlated wave function and optional
         !! otherwise: the reference built here from the orbitals is a
         !! determinant expression, so it is the energy of the wave function
         !! only when the wave function is a determinant. Passing the CI energy
         !! is what keeps the correlated case checked rather than trusted.

      real(dp), allocatable :: s_mbs(:, :), mixed(:, :), projection(:, :)
      real(dp), allocatable :: valence_internal(:, :), kinetic(:, :)
      real(dp), allocatable :: kbo(:, :), interference(:, :), populations(:)
      real(dp), allocatable :: valence_density(:, :)
      real(dp), allocatable :: kin_intra(:), kin_inter(:, :)
      real(dp), allocatable :: v_atom_ao(:, :, :), v_atom_quao(:, :, :)
      real(dp), allocatable :: nuc_intra(:), nuc_inter(:, :)
      real(dp), allocatable :: full_interference(:, :), h_core(:, :), density_ao(:, :)
      type(quao_result_t) :: core_quao, full_quao
      real(dp) :: one_electron, reference, two_electron, total_energy
      real(dp), allocatable :: eri_ao(:, :, :, :), eri_quao(:, :, :, :)
      real(dp), allocatable :: two_intra(:), two_inter(:, :)
      real(dp), allocatable :: rep_inter(:, :), tot_intra(:), tot_inter(:, :)
      real(dp), allocatable :: s_ao(:, :), u_active(:, :)
      real(dp), allocatable :: cumulant(:, :, :, :), cumulant_quao(:, :, :, :)
      logical :: correlated
      real(dp) :: span_deficit, formation
      real(dp), allocatable :: free_energy(:), adaptation(:)
      real(dp), allocatable :: nuc_coulomb(:, :), two_coulomb(:, :), classical(:, :)
      integer, allocatable :: core_off(:), core_n(:), val_off(:), val_n(:)
      character(len=160) :: line
      character(len=8) :: label
      integer :: natm, iatom, i, core, valence
      logical :: loud

      if (error%has_error()) return
      loud = .true.
      if (present(verbose)) loud = verbose
      natm = size(atomic_numbers)

      ! The dimensions first: they are pure counting and they refuse the cases
      ! this cannot handle -- an element past xenon, an odd electron count --
      ! before any integral is built.
      call aambs_dimensions(atomic_numbers, n_electrons, dims, error)
      if (error%has_error()) return

      call build_aambs_molecule(atomic_numbers, element_symbols, coordinates, aambs, error)
      if (error%has_error()) return
      call aambs%overlap(s_mbs)
      call mixed_basis_overlap(mol, aambs, mixed, error)
      if (error%has_error()) return

      call mo_aambs_overlap(orbitals, mixed, s_mbs, projection, error)
      call valence_virtual_orbitals(orbitals, projection, dims, vvo, error)
      call aambs_atom_ranges(atomic_numbers, aambs, core_off, core_n, val_off, &
                             val_n, error)
      if (error%has_error()) return

      ! The valence-internal space: occupied orbitals that are not core,
      ! followed by the valence-virtual ones.
      allocate (valence_internal(size(orbitals, 1), dims%n_valocc + vvo%n_vvo))
      valence_internal(:, 1:dims%n_valocc) = &
         orbitals(:, dims%n_core + 1:dims%n_occupied)
      valence_internal(:, dims%n_valocc + 1:) = vvo%orbitals

      ! The density in the valence-internal basis. For a reference determinant
      ! this is two on the occupied diagonal and zero elsewhere, which is what
      ! `quasi_atomic_orbitals` assumes when nothing is passed. For a correlated
      ! one it is neither diagonal nor idempotent, and it has to be projected:
      !
      !     D_vi = W^T diag(n) W ,   W = C^T S V
      !
      ! with `V` the valence-internal orbitals and `C` the natural orbitals.
      ! **The projection is not optional even though the density is diagonal in
      ! `C`.** The valence-virtual orbitals are combinations of virtuals with
      ! *different* occupations, so the density restricted to this space has
      ! off-diagonal elements even when the full one does not. Reading
      ! occupations straight off the diagonal would silently drop them.
      if (present(occupations)) then
         call project_density(mol, orbitals, occupations, valence_internal, &
                              valence_density, error)
         if (error%has_error()) return
         call quasi_atomic_orbitals(atomic_numbers, valence_internal, dims%n_valocc, &
                                    mixed, val_off, val_n, quao, error, &
                                    valence_density=valence_density)
      else
         call quasi_atomic_orbitals(atomic_numbers, valence_internal, dims%n_valocc, &
                                    mixed, val_off, val_n, quao, error)
      end if
      call orient_quasi_atomic_orbitals(quao, error)
      if (error%has_error()) return

      call mol%kinetic(kinetic)
      call kinetic_bond_orders(quao, kinetic, kbo, error, interference)
      call label_quasi_atomic_orbitals(quao, mol, aambs, mixed, s_mbs, &
                                       atomic_numbers, coordinates, labels, error)
      if (error%has_error()) return

      if (loud) then
         ! Populations first: the headline of the analysis, and the number a
         ! reader checks against chemical intuition before looking at anything
         ! else. The core electrons are added back because the construction only
         ! ever sees the valence.
         allocate (populations(natm))
         populations = 0.0_dp
         do i = 1, quao%n_quao
            populations(quao%atom_of(i)) = populations(quao%atom_of(i)) &
                                           + quao%population_bond_order(i, i)
         end do
         do iatom = 1, natm
            call aambs_element_counts(atomic_numbers(iatom), core, valence, error)
            if (error%has_error()) return
            populations(iatom) = populations(iatom) + 2.0_dp*real(core, dp)
         end do

         call logger%info("")
         call logger%info("  ==== quasi-atomic bonding analysis "// &
                          "================================")
         call logger%info("")
         call logger%info("  atoms")
         call logger%info("     atom      population      charge")
         do iatom = 1, natm
            ! Into a fixed-width buffer rather than straight through an `i0`:
            ! the atom number is variable width, so writing it inline shifts
            ! every column right the moment a molecule has ten atoms.
            write (label, "(a,i0)") trim(adjustl(element_symbols(iatom)))//" ", iatom
            write (line, "(4x,a8,f13.6,f12.6)") label, populations(iatom), &
               real(atomic_numbers(iatom), dp) - populations(iatom)
            call logger%info(trim(line))
         end do
         write (line, "(4x,a8,f13.6)") "total   ", sum(populations)
         call logger%info(trim(line))

         ! **The populations need not sum to the electron count, and when they
         ! do not it is a statement about the analysis rather than an error.**
         ! The quasi-atomic orbitals span the occupied valence orbitals plus the
         ! valence-virtual ones, and nothing else. A reference determinant puts
         ! every valence electron inside that space, so the sum is exact. A
         ! correlated density does not: the valence-virtual orbitals are chosen
         ! to look like free-atom orbitals, not to be natural orbitals, so
         ! whatever occupation sits outside their span is not recovered.
         !
         ! It is small -- 3.7 millielectrons of 14 for N2 CAS(6,6) in cc-pVDZ --
         ! and it is the honest measure of how much of the correlated wave
         ! function this analysis is describing. Printed rather than hidden,
         ! because a population analysis quietly losing electrons is the kind of
         ! thing that gets discovered much later and from much further away.
         if (abs(sum(populations) - real(n_electrons, dp)) > 1.0e-6_dp) then
            write (line, "(a,f10.6,a,i0,a)") &
               "    outside the valence space ", &
               real(n_electrons, dp) - sum(populations), " of ", n_electrons, &
               " electrons"
            call logger%info(trim(line))
         end if

         call print_quao_report(.true., quao, labels, interference, element_symbols, &
                                dims%n_core, kinetic_bond_order=kbo, &
                                threshold=threshold)

         ! The same interference matrix regrouped by atom rather than by orbital
         ! pair. Unscaled, so unlike the kinetic bond orders above these are
         ! energies and they add up to one -- the valence kinetic energy, since
         ! the core is not in this basis.
         ! The decomposition needs the core, which the bonding analysis above
         ! does not: it carries most of the kinetic energy and most of the
         ! nuclear attraction, and without it the totals below would sum to a
         ! number nothing else in the calculation knows. Built by the same
         ! construction, handed the core orbitals and the core minimal-basis
         ! ranges. Left unoriented deliberately -- see `quao_interference`.
         if (dims%n_core > 0) then
            call quasi_atomic_orbitals(atomic_numbers, orbitals(:, 1:dims%n_core), &
                                       dims%n_core, mixed, core_off, core_n, &
                                       core_quao, error)
            if (error%has_error()) return
            call combine_quao_sets(core_quao, quao, full_quao, error)
         else
            full_quao = quao
         end if
         if (error%has_error()) return

         call quao_interference(full_quao, kinetic, full_interference, error)
         if (error%has_error()) return
         call kinetic_decomposition(full_quao, full_interference, natm, kin_intra, &
                                    kin_inter, error)
         if (error%has_error()) return
         call print_kinetic_decomposition(kin_intra, kin_inter, element_symbols)

         ! The other half of the one-electron energy. Split by which nucleus is
         ! attracting, so a term carries three atomic labels rather than two and
         ! an atom's density sitting in a foreign nuclear field becomes an
         ! interaction between that pair rather than a property of one atom.
         call nuclear_attraction_per_atom(mol, atomic_numbers, coordinates, v_atom_ao, error)
         if (error%has_error()) return
         call quao_nuclear_attraction(full_quao, v_atom_ao, v_atom_quao, error)
         if (error%has_error()) return
         call nuclear_decomposition(full_quao, full_quao%population_bond_order, &
                                    v_atom_quao, natm, nuc_intra, nuc_inter, error, &
                                    coulomb=nuc_coulomb)
         if (error%has_error()) return
         call print_nuclear_decomposition(nuc_intra, nuc_inter, element_symbols)

         ! **The check the decomposition exists to pass.** Everything above is
         ! internally consistent by construction -- each routine verifies its own
         ! bins add up -- but that says nothing about whether the basis they were
         ! summed in holds the whole density. This compares the one-electron
         ! energy the decomposition arrives at against the one the converged
         ! orbitals give directly, in the atomic-orbital basis, with no
         ! quasi-atomic anything involved.
         !
         ! For a reference determinant the two are equal to rounding: the core
         ! and valence-internal spaces together contain every occupied orbital.
         ! For a correlated density they are not, and the difference is the same
         ! occupation living outside the valence-virtual span that the population
         ! sum reports above -- a statement about the analysis rather than an
         ! error, so it is printed rather than raised.
         ! The two-electron energy, where a term carries four atomic labels
         ! rather than three and both a bra and a ket distribution can be spread
         ! across two atoms at once.
         call mol%eris(eri_ao)
         call quao_eris(full_quao, eri_ao, eri_quao, error)
         if (error%has_error()) return

         ! What correlation adds. The determinant expression below is evaluated
         ! over the whole quasi-atomic basis whatever the wave function; the
         ! cumulant is the part no determinant can produce, it is zero unless
         ! all four of its indices are active, and it is carried here from the
         ! orbitals it was computed in rather than from the natural ones.
         correlated = present(active_dm2) .and. present(active_dm1) .and. &
                      present(active_orbitals)
         if (correlated) then
            call mol%overlap(s_ao)
            call quao_projection(full_quao, s_ao, active_orbitals, u_active, &
                                 span_deficit, error)
            if (error%has_error()) return
            ! The active orbitals are not quite inside the span -- they were
            ! optimised to lower an energy, the valence-virtual ones to look
            ! atomic. Reported for the same reason the population shortfall is:
            ! it says how much of the correlation this is describing.
            if (span_deficit > 1.0e-6_dp) then
               write (line, "(a,es10.2)") &
                  "    active orbitals outside the quasi-atomic span ", span_deficit
               call logger%info(trim(line))
            end if
            call active_cumulant(active_dm1, active_dm2, cumulant)
            call transform_cumulant(u_active, cumulant, cumulant_quao, error)
            if (error%has_error()) return
            call two_electron_decomposition(full_quao, full_quao%population_bond_order, &
                                            eri_quao, natm, two_intra, two_inter, &
                                            two_electron, error, cumulant=cumulant_quao, &
                                            coulomb=two_coulomb)
         else
            call two_electron_decomposition(full_quao, full_quao%population_bond_order, &
                                            eri_quao, natm, two_intra, two_inter, &
                                            two_electron, error, coulomb=two_coulomb)
         end if
         if (error%has_error()) return
         call print_two_electron_decomposition(two_intra, two_inter, element_symbols)

         ! **The check the whole thing exists to pass.** Each routine above
         ! verifies its own bins add up, which says nothing about whether the
         ! basis they were summed in holds the whole density. This compares the
         ! energy the decomposition arrives at against the energy of the same
         ! wave function computed in the atomic-orbital basis, with nothing
         ! quasi-atomic anywhere in it.
         !
         ! For a reference determinant the two are equal to rounding: the core
         ! and valence-internal spaces together contain every occupied orbital.
         ! For a correlated density they are not, and the difference is the same
         ! occupation living outside the valence-virtual span the population sum
         ! reports above -- a statement about the analysis rather than an error,
         ! so it is printed rather than raised.
         ! The nuclear repulsion needs no decomposing -- it is a sum over pairs
         ! already -- so adding it in the same convention resolves the whole
         ! energy into atoms and pairs with nothing left outside.
         call nuclear_repulsion_pairs(atomic_numbers, coordinates, rep_inter, error)
         if (error%has_error()) return
         allocate (tot_intra(natm), tot_inter(natm, natm))
         tot_intra = kin_intra + nuc_intra + two_intra
         tot_inter = kin_inter + nuc_inter + two_inter + rep_inter
         call print_total_decomposition(tot_intra, tot_inter, element_symbols)

         ! Everything an electrostatic model could account for: each atom's
         ! density in the other's nuclear field, the repulsion of their two
         ! charge clouds, and the repulsion of their nuclei. The rest is
         ! interference, and the kinetic contribution is entirely interference
         ! by construction -- two orbitals on different atoms is what the word
         ! means.
         allocate (classical(natm, natm))
         classical = nuc_coulomb + two_coulomb + rep_inter
         call print_interatomic_split(tot_inter, classical, element_symbols)

         ! The atoms as they would be on their own. Subtracting them turns a
         ! table of large numbers into the energy of formation, which is the
         ! quantity chemistry is about: what it costs to prepare each atom in
         ! the shape the molecule needs, against what the pairs give back.
         call free_atom_energies(mol, free_energy, error)
         if (error%has_error()) return
         allocate (adaptation(natm))
         adaptation = tot_intra - free_energy
         formation = sum(adaptation) + 0.5_dp*sum(tot_inter)
         call print_formation(adaptation, free_energy, tot_intra, tot_inter, &
                              element_symbols, formation)

         ! Subtracting a constant from every atomic term cannot change what the
         ! pieces sum to, so this must equal the total less the free atoms. It
         ! catches a free-atom energy landing on the wrong atom, which nothing
         ! else here would notice.
         if (abs(formation - (kinetic_total(tot_intra, tot_inter) - sum(free_energy))) &
             > 1.0e-8_dp) then
            call error%set(ERROR_VALIDATION, "the energy of formation does not match "// &
                           "the total less the free atoms, so the atomic references "// &
                           "were not subtracted from the atoms they belong to.")
            return
         end if

         one_electron = kinetic_total(kin_intra, kin_inter) + &
                        kinetic_total(nuc_intra, nuc_inter)
         call mol%core_hamiltonian(h_core)
         call one_electron_reference(orbitals, dims%n_occupied, occupations, &
                                     h_core, density_ao, reference)
         total_energy = one_electron + two_electron + mol%nuclear_repulsion()
         ! The reference built from the orbitals is a determinant expression and
         ! describes a correlated wave function not at all, so for one of those
         ! the energy the calculation reported is what the decomposition is held
         ! against instead. Either way it is a number reached without any of
         ! this arithmetic, which is the whole point of comparing.
         if (correlated .and. .not. present(reference_energy)) then
            call error%set(ERROR_VALIDATION, "a correlated energy decomposition was "// &
                           "asked for without the energy to check it against. The "// &
                           "reference built here assumes a single determinant and "// &
                           "would disagree by the whole correlation energy.")
            return
         end if
         if (present(reference_energy)) then
            reference = reference_energy
         else
            reference = reference + two_electron_reference(density_ao, eri_ao) + &
                        mol%nuclear_repulsion()
         end if

         call logger%info("")
         write (line, "(4x,a30,f16.6,a)") "one-electron", one_electron, " hartree"
         call logger%info(trim(line))
         write (line, "(4x,a30,f16.6,a)") "two-electron", two_electron, " hartree"
         call logger%info(trim(line))
         write (line, "(4x,a30,f16.6,a)") "nuclear repulsion", &
            mol%nuclear_repulsion(), " hartree"
         call logger%info(trim(line))
         write (line, "(4x,a30,f16.6,a)") "total", total_energy, " hartree"
         call logger%info(trim(line))
         if (present(reference_energy)) then
            write (line, "(4x,a30,f16.6,a)") "reported by the calculation", reference, &
               " hartree"
         else
            write (line, "(4x,a30,f16.6,a)") "from the orbitals directly", reference, &
               " hartree"
         end if
         call logger%info(trim(line))
         if (abs(kinetic_total(tot_intra, tot_inter) - total_energy) > 1.0e-8_dp) then
            call error%set(ERROR_VALIDATION, "the atom and atom-pair table does not "// &
                           "sum to the total energy, so the four contributions were "// &
                           "not accumulated consistently.")
            return
         end if
         if (abs(total_energy - reference) > 1.0e-8_dp) then
            write (line, "(4x,a30,f16.6,a)") "outside the quasi-atomic span", &
               reference - total_energy, " hartree"
            call logger%info(trim(line))
         end if

         ! Two diagnostics, at the end rather than the top: they say whether to
         ! believe the tables above, which is not a question worth asking until
         ! the tables have been read. The atomic character is how much of each
         ! orbital stayed on its own atom -- Paper II's functional, one when the
         ! orbitals are exactly free-atom ones. The gap is the separation the
         ! valence-virtual selection cut through; Paper I reports 0.99999
         ! against 0.105-0.272, and anything narrower means the valence space
         ! was not cleanly separable and the whole analysis is on sand.
         call logger%info("")
         write (line, "(a,f8.4,a,i0,a)") "    atomic character   ", &
            quao%atomic_character, "   (", quao%sweeps, " refinement sweeps)"
         call logger%info(trim(line))
         write (line, "(a,f8.4,a,f8.4)") "    valence gap        ", &
            vvo%smallest_retained, "   against rejected ", vvo%largest_rejected
         call logger%info(trim(line))
         call logger%info("")
         deallocate (populations)
      end if

      call aambs%destroy()
      deallocate (s_mbs, mixed, projection, valence_internal, kinetic, kbo)
      deallocate (interference, core_off, core_n, val_off, val_n)
   end subroutine run_quao_analysis

   subroutine project_density(mol, natural, occupations, valence_internal, &
                              density, error)
      !! `D_vi = W^T diag(n) W` with `W = C^T S V`
      !!
      !! The one-particle density, expressed in the valence-internal basis. `C`
      !! must be natural orbitals -- the density diagonal in that basis with the
      !! given occupations -- because that is the assumption the whole
      !! expression rests on.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: natural(:, :)           !! (n_ao, n_mo)
      real(dp), intent(in) :: occupations(:)          !! (n_mo)
      real(dp), intent(in) :: valence_internal(:, :)  !! (n_ao, n_val)
      real(dp), allocatable, intent(out) :: density(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: overlap(:, :), work(:, :), weights(:, :), scaled(:, :)
      integer :: n_ao, n_mo, n_val, i

      if (error%has_error()) return
      n_ao = size(natural, 1)
      n_mo = size(natural, 2)
      n_val = size(valence_internal, 2)
      if (size(occupations) /= n_mo) then
         call error%set(ERROR_VALIDATION, "there are "//to_char(size(occupations))// &
                        " occupations for "//to_char(n_mo)//" orbitals.")
         return
      end if

      call mol%overlap(overlap)
      allocate (work(n_ao, n_val), weights(n_mo, n_val))
      call pic_gemm(overlap, valence_internal, work)
      call pic_gemm(natural, work, weights, transa="T")

      allocate (scaled(n_mo, n_val), density(n_val, n_val))
      do i = 1, n_mo
         scaled(i, :) = occupations(i)*weights(i, :)
      end do
      call pic_gemm(weights, scaled, density, transa="T")

      deallocate (overlap, work, weights, scaled)
   end subroutine project_density

   subroutine one_electron_reference(orbitals, n_occupied, occupations, h, &
                                     density_ao, energy)
      !! The one-electron energy straight from the converged orbitals
      !!
      !! Deliberately built without touching anything quasi-atomic, so that
      !! comparing it against the decomposition is a real check rather than a
      !! restatement. `sum(gamma * h)` rather than a trace with a transpose,
      !! since both matrices are symmetric.
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_occupied
      real(dp), intent(in), optional :: occupations(:)
      real(dp), intent(in) :: h(:, :)
      real(dp), allocatable, intent(out) :: density_ao(:, :)
      real(dp), intent(out) :: energy

      integer :: n_ao, i, n_use

      n_ao = size(orbitals, 1)
      allocate (density_ao(n_ao, n_ao))
      density_ao = 0.0_dp

      if (present(occupations)) then
         n_use = min(size(occupations), size(orbitals, 2))
         do i = 1, n_use
            if (abs(occupations(i)) < 1.0e-14_dp) cycle
            density_ao = density_ao + occupations(i)* &
                         spread(orbitals(:, i), 2, n_ao)*spread(orbitals(:, i), 1, n_ao)
         end do
      else
         do i = 1, n_occupied
            density_ao = density_ao + 2.0_dp* &
                         spread(orbitals(:, i), 2, n_ao)*spread(orbitals(:, i), 1, n_ao)
         end do
      end if

      energy = sum(density_ao*h)
   end subroutine one_electron_reference

   pure function two_electron_reference(density_ao, eri_ao) result(energy)
      !! The two-electron energy straight from the atomic-orbital density
      !!
      !! The same `Gamma = gamma gamma - (1/2) gamma gamma` as the decomposition
      !! uses, contracted in the basis the integrals were computed in. Built
      !! independently so that comparing the two is a check and not a
      !! restatement.
      real(dp), intent(in) :: density_ao(:, :)
      real(dp), intent(in) :: eri_ao(:, :, :, :)
      real(dp) :: energy

      integer :: n, mu, nu, la, si

      n = size(density_ao, 1)
      energy = 0.0_dp
      do mu = 1, n
         do nu = 1, n
            do la = 1, n
               do si = 1, n
                  energy = energy + 0.5_dp*(density_ao(mu, nu)*density_ao(la, si) &
                                            - 0.5_dp*density_ao(mu, si)*density_ao(la, nu)) &
                           *eri_ao(mu, nu, la, si)
               end do
            end do
         end do
      end do
   end function two_electron_reference

   subroutine print_formation(adaptation, free_energy, intra, inter, &
                              element_symbols, formation)
      !! What it costs to prepare the atoms, against what the pairs give back
      !!
      !! **The adaptation energy is positive and that is not a bug.** An atom in
      !! a molecule is deformed -- promoted, and in a polar bond stripped of
      !! charge -- and both cost energy against the free atom. Binding comes
      !! from the interatomic terms being more negative than the adaptation is
      !! positive, which is why the two are printed side by side.
      real(dp), intent(in) :: adaptation(:), free_energy(:), intra(:), inter(:, :)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: formation

      character(len=160) :: line
      character(len=16) :: label
      integer :: natm, a

      natm = size(adaptation)
      call logger%info("")
      call logger%info("  energy of formation")
      call logger%info("     atom          in molecule      free atom     adaptation")
      do a = 1, natm
         write (label, "(a,i0)") trim(adjustl(element_symbols(a)))//" ", a
         write (line, "(4x,a10,3f15.6)") label, intra(a), free_energy(a), adaptation(a)
         call logger%info(trim(line))
      end do
      call logger%info("")
      write (line, "(4x,a24,f16.6,a)") "adaptation total", sum(adaptation), " hartree"
      call logger%info(trim(line))
      write (line, "(4x,a24,f16.6,a)") "interatomic total", 0.5_dp*sum(inter), " hartree"
      call logger%info(trim(line))
      write (line, "(4x,a24,f16.6,a,f12.3,a)") "energy of formation", formation, &
         " hartree", formation*HARTREE_TO_KCALMOL, " kcal/mol"
      call logger%info(trim(line))
   end subroutine print_formation

end module mqc_libcint_bonding
