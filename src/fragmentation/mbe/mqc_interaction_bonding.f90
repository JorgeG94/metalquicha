!! The reference fragment's bonding with the rest of each term it is in
module mqc_interaction_bonding
   !! Under `driver: "InteractionEnergy"` every term holding the reference
   !! fragment is computed, and with `properties.bonding_analysis` each one
   !! returns its quasi-atomic bonding tables. This module keeps, per term, only
   !! what crosses between the reference and the other monomers -- the part
   !! that describes their interaction -- renumbers it to the system's atoms,
   !! and prints it.
   !!
   !! The tables are per term and are not combined across terms: a kinetic
   !! bond order between two atoms is a property of the calculation it came
   !! from, and there is no many-body expansion of one.
   use pic_types, only: dp, int64, int_index
   use pic_io, only: to_char
   use pic_logger, only: logger => global_logger
   use pic_sorting, only: sort_index
   use mqc_combinatorics, only: is_auxiliary_row, real_count_of
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_elements, only: element_number_to_symbol
   use mqc_quao_rows, only: quao_rows_t, quao_type_name, QUAO_ROW_BOND, QUAO_TYPE_SIGMA, &
                            QUAO_TYPE_PI
   use mqc_result_types, only: calculation_result_t
   use mqc_json_output_types, only: interaction_bonding_term_t
   use mqc_physical_fragment, only: system_geometry_t, physical_fragment_t, &
                                    build_fragment_from_indices, get_monomer_atom_list
   implicit none
   private

   public :: collect_interaction_bonding
   public :: print_interaction_bonding

contains

   subroutine collect_interaction_bonding(polymers, fragment_count, reference, results, &
                                          sys_geom, terms, error)
      !! One entry per real term holding the reference and another monomer
      !! whose result carries bonding rows, in term-list order
      !!
      !! Unallocated `terms` means no such term had rows. A term whose
      !! analysis failed is skipped; its calculation already warned.
      ! TODO(mqc): the analysis ran on every term, including the subsets that
      ! do not hold the reference, whose rows are discarded here. The worker
      ! has the monomer list but `physical_fragment_t` does not, so the backend
      ! cannot skip them; with `energy_decomposition` on that is a dense
      ! two-electron transformation per discarded term.
      integer, intent(in) :: polymers(:, :)
      integer(int64), intent(in) :: fragment_count
      integer, intent(in) :: reference
         !! Monomer number, 1-based
      type(calculation_result_t), intent(in) :: results(:)
      type(system_geometry_t), intent(in) :: sys_geom
      type(interaction_bonding_term_t), allocatable, intent(out) :: terms(:)
      type(error_t), intent(out) :: error

      integer, allocatable :: monomer_of(:), atoms(:), picked(:)
      integer :: m, n_here, n_terms, t
      integer(int64) :: i

      ! Which monomer each system atom belongs to, 1-based both ways.
      allocate (monomer_of(sys_geom%total_atoms))
      monomer_of = 0
      do m = 1, sys_geom%n_monomers
         call get_monomer_atom_list(sys_geom, m, atoms, n_here)
         if (n_here > 0) monomer_of(atoms + 1) = m
      end do

      allocate (picked(fragment_count))
      n_terms = 0
      do i = 1_int64, fragment_count
         if (is_auxiliary_row(polymers(i, :))) cycle
         if (.not. any(polymers(i, :) == reference)) cycle
         if (real_count_of(polymers(i, :)) < 2) cycle
         if (.not. results(i)%has_quao_rows) cycle
         n_terms = n_terms + 1
         picked(n_terms) = int(i)
      end do
      if (n_terms == 0) return

      allocate (terms(n_terms))
      do t = 1, n_terms
         call crossing_rows(polymers(picked(t), :), picked(t), reference, &
                            results(picked(t))%quao_rows, sys_geom, monomer_of, terms(t), error)
         if (error%has_error()) return
      end do
   end subroutine collect_interaction_bonding

   subroutine crossing_rows(row, term, reference, rows, sys_geom, monomer_of, out, error)
      !! Keep the orbital and atom pairs of one term that cross between the
      !! reference and the other monomers, renumbered to the system
      integer, intent(in) :: row(:)
         !! This term's entry in the term list, zero-padded, negative for a ghost
      integer, intent(in) :: term
      integer, intent(in) :: reference
      type(quao_rows_t), intent(in) :: rows
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: monomer_of(:)
      type(interaction_bonding_term_t), intent(out) :: out
      type(error_t), intent(out) :: error

      type(physical_fragment_t) :: fragment
      integer, allocatable :: to_system(:), owner(:), kept(:)
      logical, allocatable :: usable(:)
      real(dp), allocatable :: key(:)
      integer(int_index), allocatable :: order(:)
      integer :: n_atoms, n_real, a, b, k, n_kept, n_pairs, p, cap

      out%term = term
      out%monomers = pack(row, row /= 0)
      if (allocated(sys_geom%bonds)) then
         call build_fragment_from_indices(sys_geom, out%monomers, fragment, error, sys_geom%bonds)
      else
         call build_fragment_from_indices(sys_geom, out%monomers, fragment, error)
      end if
      if (error%has_error()) return

      ! Local atom to system atom, and to the monomer that owns it. A cap
      ! takes the owner of the atom it is bonded to, so a row reaching a cap
      ! can still be placed on one side; it is then dropped, since a cap is not
      ! an atom of the system.
      n_atoms = fragment%n_atoms
      n_real = n_atoms - fragment%n_caps
      allocate (to_system(n_atoms), owner(n_atoms), usable(n_atoms))
      to_system = 0
      to_system(1:n_real) = fragment%local_to_global(1:n_real)
      usable = .false.
      usable(1:n_real) = .true.
      if (allocated(fragment%is_ghost)) usable = usable .and. .not. fragment%is_ghost
      owner = 0
      do a = 1, n_real
         owner(a) = monomer_of(to_system(a))
      end do
      do cap = 1, fragment%n_caps
         owner(n_real + cap) = monomer_of(fragment%cap_bonded_to(cap) + 1)
      end do
      if (rows%n > 0) then
         if (maxval(rows%atom) > n_atoms) then
            call error%set(ERROR_VALIDATION, "the bonding rows of term "//to_char(term)// &
                           " name an atom past the fragment's "//to_char(n_atoms))
            return
         end if
      end if

      ! ---- orbital pairs ---------------------------------------------------
      allocate (kept(max(rows%n, 1)))
      n_kept = 0
      do k = 1, rows%n
         a = rows%atom(1, k)
         b = rows%atom(2, k)
         if ((owner(a) == reference) .eqv. (owner(b) == reference)) cycle
         if (.not. (usable(a) .and. usable(b))) then
            out%omitted_rows = out%omitted_rows + 1
            cycle
         end if
         n_kept = n_kept + 1
         kept(n_kept) = k
      end do
      call copy_rows(rows, kept(1:n_kept), to_system, owner, usable, out)

      ! ---- atom pairs --------------------------------------------------------
      ! Every reference atom against every other real atom of the term, kept
      ! at the same threshold the orbital pairs were, strongest bonding first.
      allocate (out%pair_atoms(2, n_real*n_real), out%pair_bond_index(n_real*n_real), &
                out%pair_kinetic_bond_order(n_real*n_real))
      n_pairs = 0
      if (allocated(rows%atom_kinetic_bond_order)) then
         do a = 1, n_real
            if (.not. usable(a) .or. owner(a) /= reference) cycle
            do b = 1, n_real
               if (.not. usable(b) .or. owner(b) == reference) cycle
               if (abs(rows%atom_kinetic_bond_order(a, b)) < rows%threshold) cycle
               n_pairs = n_pairs + 1
               out%pair_atoms(:, n_pairs) = [to_system(a), to_system(b)]
               out%pair_bond_index(n_pairs) = rows%atom_bond_index(a, b)
               out%pair_kinetic_bond_order(n_pairs) = rows%atom_kinetic_bond_order(a, b)
            end do
         end do
      end if
      allocate (key(n_pairs), order(n_pairs))
      key = out%pair_kinetic_bond_order(1:n_pairs)
      if (n_pairs > 1) call sort_index(key, order)
      if (n_pairs == 1) order(1) = 1
      out%pair_atoms = out%pair_atoms(:, [(int(order(p)), p=1, n_pairs)])
      out%pair_bond_index = out%pair_bond_index([(int(order(p)), p=1, n_pairs)])
      out%pair_kinetic_bond_order = out%pair_kinetic_bond_order([(int(order(p)), p=1, n_pairs)])

      call fragment%destroy()
   end subroutine crossing_rows

   subroutine copy_rows(rows, kept, to_system, owner, usable, out)
      !! The kept orbital pairs, renumbered to system atoms
      type(quao_rows_t), intent(in) :: rows
      integer, intent(in) :: kept(:)
      integer, intent(in) :: to_system(:), owner(:)
      logical, intent(in) :: usable(:)
      type(interaction_bonding_term_t), intent(inout) :: out

      integer :: n, k, side, partner

      n = size(kept)
      out%rows%n = n
      out%rows%threshold = rows%threshold
      out%rows%orientation_stalled = rows%orientation_stalled
      out%rows%kind = rows%kind(kept)
      out%rows%orbital = rows%orbital(:, kept)
      out%rows%orbital_type = rows%orbital_type(:, kept)
      out%rows%dominant_l = rows%dominant_l(:, kept)
      out%rows%occupation = rows%occupation(:, kept)
      out%rows%bond_order = rows%bond_order(kept)
      out%rows%kinetic_bond_order = rows%kinetic_bond_order(kept)
      allocate (out%rows%atom(2, n), out%rows%partner_atom(2, n), out%monomer_of_atom(2, n))
      do k = 1, n
         do side = 1, 2
            out%rows%atom(side, k) = to_system(rows%atom(side, kept(k)))
            out%monomer_of_atom(side, k) = owner(rows%atom(side, kept(k)))
            partner = rows%partner_atom(side, kept(k))
            out%rows%partner_atom(side, k) = 0
            if (partner > 0) then
               if (usable(partner)) out%rows%partner_atom(side, k) = to_system(partner)
            end if
         end do
      end do
   end subroutine copy_rows

   subroutine print_interaction_bonding(terms, sys_geom, reference)
      !! The crossing rows of every term, as a table per term
      !!
      !! Atom and monomer numbers are 0-based, as the deck's are.
      type(interaction_bonding_term_t), intent(in) :: terms(:)
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: reference
         !! Monomer number, 1-based

      character(len=200) :: line
      character(len=24) :: left, right
      character(len=12) :: direction
      integer :: t, k, a, b

      if (size(terms) == 0) return
      call logger%info("")
      write (line, "(a,i0,a,f0.2,a)") "Bonding of reference fragment ", reference - 1, &
         " with its environment (quasi-atomic, |kinetic bond order| >= ", &
         terms(1)%rows%threshold, " kcal/mol)"
      call logger%info(trim(line))
      call logger%info("  atoms and monomers are 0-based; reference atoms are marked *")

      do t = 1, size(terms)
         call logger%info("")
         line = "  term "//to_char(terms(t)%term)//": monomers"
         do k = 1, size(terms(t)%monomers)
            line = trim(line)//" "//to_char(abs(terms(t)%monomers(k)) - 1)
         end do
         call logger%info(trim(line))

         call logger%info("     atom pair                  index    kcal/mol")
         if (size(terms(t)%pair_kinetic_bond_order) == 0) then
            call logger%info("     (none above threshold)")
         end if
         do k = 1, size(terms(t)%pair_kinetic_bond_order)
            a = terms(t)%pair_atoms(1, k)
            b = terms(t)%pair_atoms(2, k)
            left = atom_label(sys_geom, a)//"*"
            right = atom_label(sys_geom, b)
            write (line, "(5x,a10,a3,a10,f9.4,f12.2)") left, " - ", right, &
               terms(t)%pair_bond_index(k), terms(t)%pair_kinetic_bond_order(k)
            call logger%info(trim(line))
         end do

         call logger%info("     orbital pair, donor first"//repeat(" ", 27)// &
                          "direction       order    kcal/mol")
         if (terms(t)%rows%n == 0) call logger%info("     (none above threshold)")
         do k = 1, terms(t)%rows%n
            left = orbital_label(sys_geom, terms(t), k, 1, reference)
            right = orbital_label(sys_geom, terms(t), k, 2, reference)
            if (terms(t)%rows%kind(k) == QUAO_ROW_BOND) then
               direction = "bond"
            else if (terms(t)%monomer_of_atom(1, k) == reference) then
               direction = "ref -> env"
            else
               direction = "env -> ref"
            end if
            write (line, "(5x,a24,a4,a24,a12,f9.4,f12.2)") left, " -- ", right, &
               direction, terms(t)%rows%bond_order(k), terms(t)%rows%kinetic_bond_order(k)
            call logger%info(trim(line))
         end do
         if (terms(t)%omitted_rows > 0) then
            call logger%info("     ("//to_char(terms(t)%omitted_rows)// &
                             " orbital pairs reaching a cap or ghost atom not shown)")
         end if
         if (terms(t)%rows%orientation_stalled) then
            call logger%info("     (orientation stopped at its sweep limit; atom pairs are "// &
                             "unaffected, orbital pairs resolved to its final angle)")
         end if
      end do
      call logger%info("")
   end subroutine print_interaction_bonding

   function atom_label(sys_geom, atom) result(label)
      !! "O 30" for system atom `atom`, 1-based in and 0-based out
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: atom
      character(len=:), allocatable :: label

      label = trim(element_number_to_symbol(sys_geom%element_numbers(atom)))//" "// &
              to_char(atom - 1)
   end function atom_label

   function orbital_label(sys_geom, term, k, side, reference) result(label)
      !! "O 30* p-lone" or "N 4-H 5 sigma" for one end of one orbital pair
      type(system_geometry_t), intent(in) :: sys_geom
      type(interaction_bonding_term_t), intent(in) :: term
      integer, intent(in) :: k, side, reference
      character(len=:), allocatable :: label

      character(len=10) :: kind
      integer :: orbital_type

      label = atom_label(sys_geom, term%rows%atom(side, k))
      if (term%monomer_of_atom(side, k) == reference) label = label//"*"
      ! A bonding orbital is named by its bond, as the per-calculation report
      ! names it; anything else by its atom alone.
      orbital_type = term%rows%orbital_type(side, k)
      if ((orbital_type == QUAO_TYPE_SIGMA .or. orbital_type == QUAO_TYPE_PI) .and. &
          term%rows%partner_atom(side, k) > 0) then
         label = label//"-"//atom_label(sys_geom, term%rows%partner_atom(side, k))
      end if
      kind = quao_type_name(orbital_type, term%rows%dominant_l(side, k))
      label = label//" "//trim(kind)
   end function orbital_label

end module mqc_interaction_bonding
