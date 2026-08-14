!! The fragment molecular orbital method, two-body, for non-covalent fragments
module mqc_libcint_fmo
   !! FMO2: every fragment solved in the electrostatic field of all the others,
   !! then every pair solved the same way, with the pair energies correcting the
   !! sum of the monomers.
   !!
   !! **What separates this from a many-body expansion.** MBE computes each
   !! fragment in vacuum and adds corrections. FMO computes each fragment in the
   !! field of the rest and iterates that field to self-consistency, so a
   !! fragment's density is polarised by its neighbours before any correction is
   !! applied. For anything with hydrogen bonds that matters a great deal: the
   !! monomers of a water cluster are substantially polarised, and an expansion
   !! that starts from unpolarised monomers has to recover all of it through
   !! its correction terms.
   !!
   !! **Non-covalent fragments only, deliberately.** Every fragment here is a
   !! whole molecule, so no bond is cut, and none of the machinery that exists
   !! to cut one -- adjusted fragment orbitals, hybrid orbital projection,
   !! hydrogen caps -- is present or needed. A partition that would sever a
   !! covalent bond is refused rather than approximated.
   !!
   !! **The energy.** With `E'` for an internal energy, meaning the fragment's
   !! own energy with its polarised density and *not* counting its interaction
   !! with the external field:
   !!
   !!     E = sum_I E'_I + sum_{I<J} (E'_IJ - E'_I - E'_J)
   !!                    + sum_{I<J} Tr(dD_IJ V_IJ)
   !!
   !! The pair term carries the whole I-J interaction, because the dimer is
   !! solved with both fragments present. The last term is the response of the
   !! pair's density to the field of everything outside it: `dD_IJ` is the
   !! dimer density minus the two monomer densities laid side by side, and it is
   !! zero exactly when the pair does not polarise.
   !!
   !! **Why two fragments is the test that matters.** With two fragments there
   !! is nothing outside the dimer, so `V_12` vanishes, the monomer terms cancel
   !! algebraically, and `E = E'_12` -- an ordinary supermolecular SCF. FMO2 on
   !! two fragments is therefore not an approximation at all, and any
   !! disagreement with a plain RHF on the same system is a bug in how fragments
   !! are built rather than an error in the method. See `validation/check_fmo`.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
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
      character(len=16) :: embedding = "mulliken"
         !! Which charges represent a fragment to its neighbours: "mulliken",
         !! "chelpg", or "none" for an unembedded expansion.
         !!
         !! Mulliken is the default because it is what FMO's point-charge
         !! approximation was defined with, but it is not the more accurate of
         !! the two at ordinary separations. Measured against supermolecular
         !! RHF on stacked water clusters of 3, 4 and 5 monomers
         !! (`validation/sweep_fmo`), CHELPG is closer in 9 of 12 cases, and
         !! every one of the three it loses is at an O-O separation of 2.70 A:
         !!
         !!     O-O sep     mulliken err    chelpg err     (5 waters, Hartree)
         !!     2.70 A         7.1e-05       -2.1e-04
         !!     2.90 A         2.0e-04       -3.1e-05
         !!     3.20 A         1.6e-04        1.3e-05
         !!     4.00 A         2.6e-05        1.1e-06
         !!
         !! That crossover is where the physics says it should be. CHELPG
         !! charges are fitted to reproduce the potential *outside* the van der
         !! Waals surface and carry no information inside it, so a neighbour
         !! pressed against that surface is being represented by charges that
         !! were never fitted for where it sits. Mulliken has no excluded
         !! region and degrades more gently there. Above about 2.9 A -- which is
         !! where hydrogen-bonded monomers actually sit -- CHELPG is better by
         !! roughly an order of magnitude, and the gap widens with separation
         !! and with cluster size.
         !!
         !! Either is far better than none: unembedded is 1 to 2 orders of
         !! magnitude worse in every case measured. And the choice is close to
         !! free, since both need the same SCF and CHELPG adds about ten percent
         !! on top of it.
      integer :: max_scc = 30
         !! Cap on the monomer self-consistency loop
      real(dp) :: scc_tol = 1.0e-6_dp
         !! Convergence on the largest change in any atomic charge
      integer :: scf_max_iter = 100
      real(dp) :: scf_energy_tol = 1.0e-9_dp
      real(dp) :: scf_density_tol = 1.0e-7_dp
   end type fmo_options_t

   type :: fmo_result_t
      !! The energy, and enough of the parts to see where it came from
      real(dp) :: energy = 0.0_dp                !! The FMO2 total
      real(dp) :: monomer_sum = 0.0_dp           !! sum_I E'_I
      real(dp) :: pair_sum = 0.0_dp              !! sum of the pair corrections
      real(dp) :: response_sum = 0.0_dp          !! sum of Tr(dD V), the last term
      integer :: scc_iterations = 0
      logical :: converged = .false.
         !! The monomer loop reached `scc_tol`, and every fragment SCF converged
      real(dp), allocatable :: monomer_energy(:)     !! E'_I
      real(dp), allocatable :: pair_correction(:, :)  !! the (I,J) term, upper triangle
      real(dp), allocatable :: charges(:)            !! final embedding charges, per atom
   end type fmo_result_t

contains

   subroutine run_fmo2(atomic_numbers, symbols, coordinates, owner, opts, res, error)
      !! Run FMO2 over a system already partitioned into whole molecules
      !!
      !! `owner(i)` is the fragment index of atom `i`, numbered 1..n_frag with
      !! no gaps -- which is what `connected_components` in
      !! [[mqc_bond_perception]] produces. Coordinates are Bohr.
      integer, intent(in) :: atomic_numbers(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: coordinates(:, :)
      integer, intent(in) :: owner(:)
      type(fmo_options_t), intent(in) :: opts
      type(fmo_result_t), intent(out) :: res
      type(error_t), intent(inout) :: error

      integer, allocatable :: members(:, :), frag_size(:), pair(:)
      real(dp), allocatable :: q(:), q_new(:), d_mon(:, :, :)
      integer, allocatable :: nao_mon(:)
      real(dp) :: e_internal, shift, e_pair, e_resp
      integer :: n_atoms, n_frag, i, j, it, nao_ij
      logical :: all_converged

      n_atoms = size(atomic_numbers)
      if (size(owner) /= n_atoms .or. size(coordinates, 2) /= n_atoms) then
         call error%set(ERROR_VALIDATION, "fmo: owner and coordinates must cover "// &
                        "every atom")
         return
      end if

      call collect_fragments(owner, n_frag, frag_size, members, error)
      if (error%has_error()) return
      if (n_frag < 2) then
         call error%set(ERROR_VALIDATION, "fmo: the system is one fragment, so there "// &
                        "is nothing to expand -- run an ordinary SCF")
         return
      end if

      allocate (res%monomer_energy(n_frag), source=0.0_dp)
      allocate (res%pair_correction(n_frag, n_frag), source=0.0_dp)
      allocate (q(n_atoms), source=0.0_dp)
      allocate (q_new(n_atoms), source=0.0_dp)
      allocate (nao_mon(n_frag), source=0)

      all_converged = .true.

      ! -- the monomer loop -------------------------------------------------
      !
      ! Each fragment sees the others as point charges, and those charges come
      ! from the previous pass. The first pass therefore has no field at all
      ! and is a set of isolated fragments; convergence is when no charge moves.
      res%converged = .false.
      do it = 1, opts%max_scc
         do i = 1, n_frag
            call fragment_scf(atomic_numbers, symbols, coordinates, &
                              members(1:frag_size(i), i), q, opts, &
                              e_internal, q_new, nao_mon(i), d_mon, i, n_frag, &
                              all_converged, error)
            if (error%has_error()) return
            res%monomer_energy(i) = e_internal
         end do

         shift = maxval(abs(q_new - q))
         q = q_new
         res%scc_iterations = it
         call logger%verbose("  fmo scc "//to_char(it)//": largest charge shift "// &
                             to_char(shift))
         if (opts%embedding == "none") exit
         if (shift < opts%scc_tol) then
            res%converged = .true.
            exit
         end if
      end do
      if (opts%embedding == "none") res%converged = .true.

      if (.not. res%converged) then
         call error%set(ERROR_VALIDATION, "fmo: the monomer charges did not settle in "// &
                        to_char(opts%max_scc)//" passes")
         return
      end if

      ! One more pass at the converged charges, to leave every monomer density
      ! consistent with the field the pairs will be computed in. Without it the
      ! densities on hand are one pass stale, and dD_IJ would carry that
      ! staleness into the response term.
      do i = 1, n_frag
         call fragment_scf(atomic_numbers, symbols, coordinates, &
                           members(1:frag_size(i), i), q, opts, &
                           e_internal, q_new, nao_mon(i), d_mon, i, n_frag, &
                           all_converged, error)
         if (error%has_error()) return
         res%monomer_energy(i) = e_internal
      end do

      res%monomer_sum = sum(res%monomer_energy)
      allocate (res%charges(n_atoms), source=q)

      ! -- the pairs ---------------------------------------------------------
      do i = 1, n_frag
         do j = i + 1, n_frag
            allocate (pair(frag_size(i) + frag_size(j)))
            pair(1:frag_size(i)) = members(1:frag_size(i), i)
            pair(frag_size(i) + 1:) = members(1:frag_size(j), j)

            call pair_scf(atomic_numbers, symbols, coordinates, pair, q, opts, &
                          d_mon(:, :, i), d_mon(:, :, j), nao_mon(i), nao_mon(j), &
                          e_internal, e_resp, nao_ij, all_converged, error)
            deallocate (pair)
            if (error%has_error()) return

            e_pair = e_internal - res%monomer_energy(i) - res%monomer_energy(j)
            res%pair_correction(i, j) = e_pair + e_resp
            res%pair_sum = res%pair_sum + e_pair
            res%response_sum = res%response_sum + e_resp
         end do
      end do

      res%energy = res%monomer_sum + res%pair_sum + res%response_sum
      res%converged = res%converged .and. all_converged
      if (.not. all_converged) then
         call error%set(ERROR_VALIDATION, "fmo: at least one fragment SCF did not "// &
                        "converge, so the total is not trustworthy")
         return
      end if
   end subroutine run_fmo2

   subroutine collect_fragments(owner, n_frag, frag_size, members, error)
      !! Turn a per-atom fragment label into a list of atoms per fragment
      integer, intent(in) :: owner(:)
      integer, intent(out) :: n_frag
      integer, allocatable, intent(out) :: frag_size(:), members(:, :)
      type(error_t), intent(inout) :: error

      integer :: i, f, widest

      n_frag = maxval(owner)
      if (minval(owner) < 1) then
         call error%set(ERROR_VALIDATION, "fmo: fragment labels start at 1")
         return
      end if

      allocate (frag_size(n_frag), source=0)
      do i = 1, size(owner)
         frag_size(owner(i)) = frag_size(owner(i)) + 1
      end do
      if (any(frag_size == 0)) then
         call error%set(ERROR_VALIDATION, "fmo: fragment labels have a gap, so one "// &
                        "fragment has no atoms")
         return
      end if

      widest = maxval(frag_size)
      allocate (members(widest, n_frag), source=0)
      frag_size = 0
      do i = 1, size(owner)
         f = owner(i)
         frag_size(f) = frag_size(f) + 1
         members(frag_size(f), f) = i
      end do
   end subroutine collect_fragments

   subroutine fragment_scf(z, symbols, coords, atoms, q, opts, e_internal, q_all, &
                           nao, d_store, which, n_frag, all_converged, error)
      !! One fragment in the field of the charges on every atom outside it
      !!
      !! `q_all` is the whole system's charge vector and only this fragment's
      !! entries are written. It is passed whole rather than as a section
      !! because a vector-subscripted section cannot be an output argument.
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: atoms(:)
      real(dp), intent(in) :: q(:)
      type(fmo_options_t), intent(in) :: opts
      real(dp), intent(out) :: e_internal
      real(dp), intent(inout) :: q_all(:)
      integer, intent(inout) :: nao
      real(dp), allocatable, intent(inout) :: d_store(:, :, :)
      integer, intent(in) :: which, n_frag
      logical, intent(inout) :: all_converged
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: v(:, :), q_frag(:)

      call embedded_scf(z, symbols, coords, atoms, q, opts, mol, scf, v, error)
      if (error%has_error()) return
      if (.not. scf%converged) all_converged = .false.

      ! The internal energy: what the SCF reported, less its interaction with
      ! the external field. `h_extra` enters H linearly, so that interaction is
      ! exactly Tr(D V) and nothing else has to be unpicked.
      e_internal = scf%energy
      if (allocated(v)) e_internal = e_internal - sum(scf%density*v)

      nao = mol%nao
      allocate (q_frag(size(atoms)))
      call fragment_charges(mol, scf%density, opts, q_frag, error)
      if (error%has_error()) return
      q_all(atoms) = q_frag

      ! The densities are kept for the response term, in a (nao_max, nao_max,
      ! n_frag) box because the fragments need not be the same size.
      if (.not. allocated(d_store)) then
         allocate (d_store(mol%nao, mol%nao, n_frag), source=0.0_dp)
      else if (size(d_store, 1) < mol%nao) then
         call regrow(d_store, mol%nao, n_frag)
      end if
      d_store(1:mol%nao, 1:mol%nao, which) = scf%density
   end subroutine fragment_scf

   subroutine pair_scf(z, symbols, coords, atoms, q, opts, d_i, d_j, nao_i, nao_j, &
                       e_internal, e_response, nao, all_converged, error)
      !! One pair in the field of everything outside the pair
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: atoms(:)
      real(dp), intent(in) :: q(:)
      type(fmo_options_t), intent(in) :: opts
      real(dp), intent(in) :: d_i(:, :), d_j(:, :)
      integer, intent(in) :: nao_i, nao_j
      real(dp), intent(out) :: e_internal, e_response
      integer, intent(out) :: nao
      logical, intent(inout) :: all_converged
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: v(:, :), delta(:, :)

      call embedded_scf(z, symbols, coords, atoms, q, opts, mol, scf, v, error)
      if (error%has_error()) return
      if (.not. scf%converged) all_converged = .false.

      nao = mol%nao
      e_internal = scf%energy
      e_response = 0.0_dp
      if (.not. allocated(v)) return

      e_internal = e_internal - sum(scf%density*v)

      ! The dimer is built from fragment I's atoms followed by fragment J's, and
      ! libcint orders basis functions by atom, so the monomer blocks sit
      ! contiguously at the front and back. That is what makes laying the two
      ! monomer densities into the dimer basis a pair of array slices rather
      ! than an index map -- and it is worth checking rather than assuming,
      ! because a silent mismatch here would be a plausible-looking wrong answer.
      if (nao_i + nao_j /= mol%nao) then
         call error%set(ERROR_VALIDATION, "fmo: the dimer basis is not the two monomer "// &
                        "bases end to end ("//to_char(nao_i)//" + "//to_char(nao_j)// &
                        " /= "//to_char(mol%nao)//")")
         return
      end if

      allocate (delta(mol%nao, mol%nao), source=scf%density)
      delta(1:nao_i, 1:nao_i) = delta(1:nao_i, 1:nao_i) - d_i(1:nao_i, 1:nao_i)
      delta(nao_i + 1:, nao_i + 1:) = delta(nao_i + 1:, nao_i + 1:) - d_j(1:nao_j, 1:nao_j)
      e_response = sum(delta*v)
   end subroutine pair_scf

   subroutine embedded_scf(z, symbols, coords, atoms, q, opts, mol, scf, v, error)
      !! Build a sub-molecule and solve it in the field of the outside charges
      !!
      !! `v` comes back unallocated when there is no field -- no outside atoms,
      !! or embedding turned off -- which is how the callers tell that the SCF
      !! energy needs no correction.
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: atoms(:)
      real(dp), intent(in) :: q(:)
      type(fmo_options_t), intent(in) :: opts
      type(libcint_molecule_t), intent(out) :: mol
      type(rhf_result_t), intent(out) :: scf
      real(dp), allocatable, intent(out) :: v(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: points(:, :), matrices(:, :, :), q_out(:)
      logical, allocatable :: inside(:)
      integer :: n_atoms, i, g, n_out, nelec

      n_atoms = size(z)
      call build_libcint_molecule(z(atoms), symbols(atoms), coords(:, atoms), &
                                  trim(opts%basis), mol, error)
      if (error%has_error()) return

      nelec = sum(z(atoms))
      if (mod(nelec, 2) /= 0) then
         call error%set(ERROR_VALIDATION, "fmo: fragment has an odd electron count and "// &
                        "this is a closed-shell method")
         return
      end if

      ! Everything not in this fragment, as a point charge where it is.
      allocate (inside(n_atoms), source=.false.)
      inside(atoms) = .true.
      n_out = count(.not. inside)

      if (n_out > 0 .and. opts%embedding /= "none" .and. any(abs(q) > 0.0_dp)) then
         allocate (points(3, n_out), q_out(n_out))
         g = 0
         do i = 1, n_atoms
            if (inside(i)) cycle
            g = g + 1
            points(:, g) = coords(:, i)
            q_out(g) = q(i)
         end do

         call esp_matrices(mol, points, matrices, error)
         if (error%has_error()) return

         ! An electron carries charge -1, so a positive point charge lowers its
         ! energy: the operator is -sum_g q_g/|r - R_g|.
         allocate (v(mol%nao, mol%nao), source=0.0_dp)
         do g = 1, n_out
            v = v - q_out(g)*matrices(:, :, g)
         end do
      end if

      if (allocated(v)) then
         call run_libcint_rhf(mol, nelec, opts%scf_max_iter, opts%scf_energy_tol, &
                              opts%scf_density_tol, .false., scf, error, h_extra=v)
      else
         call run_libcint_rhf(mol, nelec, opts%scf_max_iter, opts%scf_energy_tol, &
                              opts%scf_density_tol, .false., scf, error)
      end if
   end subroutine embedded_scf

   subroutine fragment_charges(mol, density, opts, q_out, error)
      !! The charges this fragment shows its neighbours
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      type(fmo_options_t), intent(in) :: opts
      real(dp), intent(out) :: q_out(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: s(:, :), q(:)

      if (opts%embedding == "none") then
         q_out = 0.0_dp
         return
      end if

      if (opts%embedding == "chelpg") then
         call chelpg_charges(mol, density, q, error)
      else
         call mol%overlap(s)
         call mulliken_charges(mol, density, s, q, error)
      end if
      if (error%has_error()) return
      q_out = q
   end subroutine fragment_charges

   subroutine regrow(box, nao, n_frag)
      !! Widen the density box when a later fragment needs more room than the
      !! first one did
      real(dp), allocatable, intent(inout) :: box(:, :, :)
      integer, intent(in) :: nao, n_frag

      real(dp), allocatable :: bigger(:, :, :)
      integer :: was

      was = size(box, 1)
      allocate (bigger(nao, nao, n_frag), source=0.0_dp)
      bigger(1:was, 1:was, :) = box
      call move_alloc(bigger, box)
   end subroutine regrow

end module mqc_libcint_fmo
