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
   use mqc_libcint_ieda, only: kinetic_decomposition, print_kinetic_decomposition, &
                               nuclear_attraction_per_atom, quao_nuclear_attraction, &
                               nuclear_decomposition, print_nuclear_decomposition
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
                                occupations)
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
      real(dp), allocatable :: s_mbs(:, :), mixed(:, :), projection(:, :)
      real(dp), allocatable :: valence_internal(:, :), kinetic(:, :)
      real(dp), allocatable :: kbo(:, :), interference(:, :), populations(:)
      real(dp), allocatable :: valence_density(:, :)
      real(dp), allocatable :: kin_intra(:), kin_inter(:, :)
      real(dp), allocatable :: v_atom_ao(:, :, :), v_atom_quao(:, :, :)
      real(dp), allocatable :: nuc_intra(:), nuc_inter(:, :)
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
         call kinetic_decomposition(quao, interference, natm, kin_intra, kin_inter, error)
         if (error%has_error()) return
         call print_kinetic_decomposition(kin_intra, kin_inter, element_symbols)

         ! The other half of the one-electron energy. Split by which nucleus is
         ! attracting, so a term carries three atomic labels rather than two and
         ! an atom's density sitting in a foreign nuclear field becomes an
         ! interaction between that pair rather than a property of one atom.
         call nuclear_attraction_per_atom(mol, atomic_numbers, coordinates, v_atom_ao, error)
         if (error%has_error()) return
         call quao_nuclear_attraction(quao, v_atom_ao, v_atom_quao, error)
         if (error%has_error()) return
         call nuclear_decomposition(quao, quao%population_bond_order, v_atom_quao, &
                                    natm, nuc_intra, nuc_inter, error)
         if (error%has_error()) return
         call print_nuclear_decomposition(nuc_intra, nuc_inter, element_symbols)

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

end module mqc_libcint_bonding
