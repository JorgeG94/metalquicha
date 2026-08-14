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
   !! *distant* fragments; applying it to all of them instead gives
   !! electrostatically embedded MBE -- a good method, and a different one.
   !! `esp = "none"` drops the embedding and leaves a plain many-body expansion.
   !! All three are here because the differences are worth measuring; see
   !! `validation/sweep_fmo`.
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
   !! **Non-covalent fragments only, deliberately.** Every fragment is a whole
   !! molecule, so no bond is cut and none of the machinery for cutting one --
   !! adjusted fragment orbitals, hybrid orbital projection, hydrogen caps -- is
   !! present or needed.
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
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_direct, only: schwarz_bounds, build_fock_direct, direct_stats_t
   use mqc_libcint_esp, only: esp_matrices
   use mqc_libcint_charges, only: ao_to_atom, mulliken_charges, chelpg_charges
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
      character(len=16) :: ptc_charges = "mulliken"
         !! Read only when `esp = "ptc"`: "mulliken" or "chelpg".
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
      type(libcint_molecule_t) :: mol
      real(dp), allocatable :: bounds(:, :)      !! Schwarz, for its own Coulomb build
      real(dp), allocatable :: density(:, :)
      integer, allocatable :: ao(:)              !! its AOs, as supersystem indices
      integer, allocatable :: atoms(:)           !! its atoms, as system indices
      integer, allocatable :: z(:)
      character(len=2), allocatable :: sym(:)
      real(dp), allocatable :: xyz(:, :)
      real(dp) :: energy = 0.0_dp                !! internal, E'
      integer :: nelec = 0
   end type fragment_t

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

      type(fragment_t), allocatable :: frag(:)
      type(libcint_molecule_t) :: super
      real(dp), allocatable :: super_bounds(:, :), j_total(:, :), q_all(:), u(:, :)
      logical, allocatable :: inside(:)
      real(dp) :: e_sum, e_prev
      integer :: n_atoms, n_frag, i, j, outer
      logical :: all_converged

      n_atoms = size(atomic_numbers)
      if (size(owner) /= n_atoms .or. size(coordinates, 2) /= n_atoms) then
         call error%set(ERROR_VALIDATION, "fmo: owner and coordinates must cover "// &
                        "every atom")
         return
      end if

      call build_fragments(atomic_numbers, symbols, coordinates, owner, opts, &
                           super, super_bounds, frag, n_frag, error)
      if (error%has_error()) return

      allocate (res%monomer_energy(n_frag), source=0.0_dp)
      allocate (res%pair_correction(n_frag, n_frag), source=0.0_dp)
      allocate (inside(n_atoms))
      all_converged = .true.

      ! -- isolated fragments, to start the outer loop from -------------------
      do i = 1, n_frag
         call inner_scf(frag(i), opts, error, all_converged=all_converged)
         if (error%has_error()) return
      end do

      ! -- the outer SCF ------------------------------------------------------
      !
      ! Each pass rebuilds the field from the densities the previous pass left,
      ! then re-solves every fragment in it. Converged when the monomer energies
      ! stop moving, which is when the densities and the field they make agree.
      if (opts%esp == "none") then
         res%converged = .true.
      else
         e_prev = sum(frag(:)%energy)
         do outer = 1, opts%max_outer
            call field_source(super, super_bounds, frag, n_frag, opts, j_total, &
                              q_all, n_atoms, error)
            if (error%has_error()) return

            do i = 1, n_frag
               inside = .false.
               inside(frag(i)%atoms) = .true.
               call embedding_operator(frag(i)%mol, frag(i)%bounds, frag(i)%ao, &
                                       frag(i)%density, inside, atomic_numbers, &
                                       coordinates, j_total, q_all, opts, u, error)
               if (error%has_error()) return
               call inner_scf(frag(i), opts, error, u, all_converged)
               if (error%has_error()) return
            end do

            e_sum = sum(frag(:)%energy)
            res%outer_iterations = outer
            res%outer_change = abs(e_sum - e_prev)
            call logger%verbose("  fmo outer "//to_char(outer)//": monomer sum "// &
                                to_char(e_sum)//", moved "//to_char(res%outer_change))
            if (res%outer_change < opts%outer_tol) then
               res%converged = .true.
               exit
            end if
            e_prev = e_sum
         end do

         if (.not. res%converged) then
            call error%set(ERROR_VALIDATION, "fmo: the outer SCF did not settle in "// &
                           to_char(opts%max_outer)//" passes; the monomer sum was still "// &
                           "moving by "//to_char(res%outer_change)//" Hartree")
            return
         end if
      end if

      do i = 1, n_frag
         res%monomer_energy(i) = frag(i)%energy
      end do
      res%monomer_sum = sum(res%monomer_energy)

      ! -- the pairs, in the field the converged monomers make ----------------
      call field_source(super, super_bounds, frag, n_frag, opts, j_total, q_all, &
                        n_atoms, error)
      if (error%has_error()) return

      do i = 1, n_frag
         do j = i + 1, n_frag
            call pair_term(frag, i, j, atomic_numbers, coordinates, j_total, q_all, &
                           opts, res, all_converged, error)
            if (error%has_error()) return
         end do
      end do

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

   subroutine build_fragments(z, symbols, coords, owner, opts, super, super_bounds, &
                              frag, n_frag, error)
      !! The supersystem, and each fragment with its place in it
      !!
      !! The supersystem exists to hold one Coulomb matrix over the whole basis.
      !! `frag%ao` records which of its functions belong to a fragment -- not a
      !! contiguous range in general, since nothing requires a fragment's atoms
      !! to be listed together.
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: owner(:)
      type(fmo_options_t), intent(in) :: opts
      type(libcint_molecule_t), intent(out) :: super
      real(dp), allocatable, intent(out) :: super_bounds(:, :)
      type(fragment_t), allocatable, intent(out) :: frag(:)
      integer, intent(out) :: n_frag
      type(error_t), intent(inout) :: error

      integer, allocatable :: super_ao_atom(:), count_per(:)
      integer :: n_atoms, i, f, k

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

      call build_libcint_molecule(z, symbols, coords, trim(opts%basis), super, error)
      if (error%has_error()) return
      call schwarz_bounds(super, super_bounds, error)
      if (error%has_error()) return
      call ao_to_atom(super, super_ao_atom)

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

         call build_libcint_molecule(frag(f)%z, frag(f)%sym, frag(f)%xyz, &
                                     trim(opts%basis), frag(f)%mol, error)
         if (error%has_error()) return
         call schwarz_bounds(frag(f)%mol, frag(f)%bounds, error)
         if (error%has_error()) return

         ! Its functions inside the supersystem, in the order the fragment
         ! itself has them: both are built atom by atom in increasing index.
         frag(f)%ao = pack([(k, k=1, super%nao)], owner(super_ao_atom) == f)
         if (size(frag(f)%ao) /= frag(f)%mol%nao) then
            call error%set(ERROR_VALIDATION, "fmo: fragment "//to_char(f)//" has "// &
                           to_char(frag(f)%mol%nao)//" basis functions on its own but "// &
                           to_char(size(frag(f)%ao))//" inside the supersystem")
            return
         end if
      end do
   end subroutine build_fragments

   subroutine field_source(super, super_bounds, frag, n_frag, opts, j_total, q_all, &
                           n_atoms, error)
      !! Whatever the embedding is built out of this pass
      !!
      !! For the exact ESP that is one Coulomb matrix over the whole basis,
      !! built from every fragment's density. It serves every fragment: the
      !! field on X is this matrix restricted to X, less the part X's own
      !! density contributed, which is a small build over X's own basis. The
      !! alternative is a supersystem build per fragment per pass, and the
      !! subtraction is exact because a four-index integral over X's functions
      !! does not care which molecule object it was computed through.
      !!
      !! For ESP-PTC it is instead one partial charge per atom.
      type(libcint_molecule_t), intent(in) :: super
      real(dp), intent(in) :: super_bounds(:, :)
      type(fragment_t), intent(in) :: frag(:)
      integer, intent(in) :: n_frag, n_atoms
      type(fmo_options_t), intent(in) :: opts
      real(dp), allocatable, intent(out) :: j_total(:, :), q_all(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: d_total(:, :), zero_h(:, :), q(:)
      type(direct_stats_t) :: stats
      integer :: f

      if (opts%esp == "ptc") then
         allocate (q_all(n_atoms), source=0.0_dp)
         do f = 1, n_frag
            call fragment_charges(frag(f), opts%ptc_charges, q, error)
            if (error%has_error()) return
            q_all(frag(f)%atoms) = q
         end do
         return
      end if
      if (opts%esp /= "exact") return

      allocate (d_total(super%nao, super%nao), source=0.0_dp)
      do f = 1, n_frag
         d_total(frag(f)%ao, frag(f)%ao) = frag(f)%density
      end do

      allocate (zero_h(super%nao, super%nao), source=0.0_dp)
      allocate (j_total(super%nao, super%nao))
      ! No core Hamiltonian and no exchange, so what comes back is J alone.
      call build_fock_direct(super, zero_h, d_total, super_bounds, j_total, stats, &
                             error, k_scale=0.0_dp, j_scale=1.0_dp)
   end subroutine field_source

   subroutine embedding_operator(mol, bounds, ao, own_density, inside, z, coords, &
                                 j_total, q_all, opts, u, error)
      !! The field the atoms marked `inside` sit in, over `mol`'s basis
      !!
      !! Works for a monomer and a dimer alike -- `inside` is what changes.
      !! Comes back unallocated when there is no field, which is how the caller
      !! knows the SCF energy needs no correction.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: bounds(:, :)
      integer, intent(in) :: ao(:)
      real(dp), intent(in) :: own_density(:, :)
      logical, intent(in) :: inside(:)
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords(:, :)
      real(dp), allocatable, intent(in) :: j_total(:, :), q_all(:)
      type(fmo_options_t), intent(in) :: opts
      real(dp), allocatable, intent(out) :: u(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: matrices(:, :, :), points(:, :), weight(:)
      real(dp), allocatable :: j_own(:, :), zero_h(:, :)
      type(direct_stats_t) :: stats
      integer :: n_atoms, i, g, n_out

      if (opts%esp == "none") return
      n_atoms = size(z)
      n_out = count(.not. inside)
      if (n_out == 0) return

      allocate (u(mol%nao, mol%nao), source=0.0_dp)
      allocate (points(3, n_out), weight(n_out))
      g = 0
      do i = 1, n_atoms
         if (inside(i)) cycle
         g = g + 1
         points(:, g) = coords(:, i)
         ! Exact ESP: the bare nuclei here, their electrons through J below.
         ! ESP-PTC: the whole neighbour atom as one partial charge, nuclei and
         ! electrons together, and no J term at all.
         if (opts%esp == "ptc") then
            weight(g) = q_all(i)
         else
            weight(g) = real(z(i), dp)
         end if
      end do

      ! An electron carries charge -1, so a positive charge lowers its energy:
      ! the operator is -sum_g w_g/|r - R_g|.
      call esp_matrices(mol, points, matrices, error)
      if (error%has_error()) return
      do g = 1, n_out
         u = u - weight(g)*matrices(:, :, g)
      end do
      if (opts%esp == "ptc") return

      ! The outside electrons, exactly: the whole-basis Coulomb matrix
      ! restricted here, less what this group's own density put into it.
      allocate (zero_h(mol%nao, mol%nao), source=0.0_dp)
      allocate (j_own(mol%nao, mol%nao))
      call build_fock_direct(mol, zero_h, own_density, bounds, j_own, stats, error, &
                             k_scale=0.0_dp, j_scale=1.0_dp)
      if (error%has_error()) return
      u = u + j_total(ao, ao) - j_own
   end subroutine embedding_operator

   subroutine pair_term(frag, a, b, z, coords, j_total, q_all, opts, res, &
                        all_converged, error)
      !! One dimer, in the field of every fragment outside it
      type(fragment_t), intent(in) :: frag(:)
      integer, intent(in) :: a, b
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords(:, :)
      real(dp), allocatable, intent(in) :: j_total(:, :), q_all(:)
      type(fmo_options_t), intent(in) :: opts
      type(fmo_result_t), intent(inout) :: res
      logical, intent(inout) :: all_converged
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: bounds(:, :), d_split(:, :), u(:, :)
      logical, allocatable :: inside(:)
      integer, allocatable :: ao(:)
      real(dp) :: e_internal, e_resp
      integer :: na, nb

      na = frag(a)%mol%nao
      nb = frag(b)%mol%nao

      call build_libcint_molecule([frag(a)%z, frag(b)%z], [frag(a)%sym, frag(b)%sym], &
                                  reshape([frag(a)%xyz, frag(b)%xyz], &
                                          [3, size(frag(a)%z) + size(frag(b)%z)]), &
                                  trim(opts%basis), mol, error)
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
      ao = [frag(a)%ao, frag(b)%ao]

      call embedding_operator(mol, bounds, ao, d_split, inside, z, coords, &
                              j_total, q_all, opts, u, error)
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

      res%pair_correction(a, b) = e_internal - frag(a)%energy - frag(b)%energy + e_resp
      res%pair_sum = res%pair_sum + e_internal - frag(a)%energy - frag(b)%energy
      res%response_sum = res%response_sum + e_resp
   end subroutine pair_term

   subroutine inner_scf(f, opts, error, u, all_converged)
      !! The inner SCF: this fragment's orbitals, against a fixed external field
      type(fragment_t), intent(inout) :: f
      type(fmo_options_t), intent(in) :: opts
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(in), optional :: u(:, :)
      logical, intent(inout), optional :: all_converged

      type(rhf_result_t) :: scf
      logical :: embedded

      embedded = .false.
      if (present(u)) embedded = allocated(u)

      if (embedded) then
         call run_libcint_rhf(f%mol, f%nelec, opts%scf_max_iter, opts%scf_energy_tol, &
                              opts%scf_density_tol, .false., scf, error, h_extra=u)
      else
         call run_libcint_rhf(f%mol, f%nelec, opts%scf_max_iter, opts%scf_energy_tol, &
                              opts%scf_density_tol, .false., scf, error)
      end if
      if (error%has_error()) return
      if (present(all_converged)) then
         if (.not. scf%converged) all_converged = .false.
      end if

      ! The internal energy: what the SCF reported, less its interaction with
      ! the field. `h_extra` enters H linearly, so that interaction is exactly
      ! Tr(D u) and nothing else has to be unpicked.
      f%energy = scf%energy
      if (embedded) f%energy = f%energy - sum(scf%density*u)
      f%density = scf%density
   end subroutine inner_scf

   subroutine fragment_charges(f, scheme, q, error)
      !! Atomic charges for a fragment whose density is already converged
      type(fragment_t), intent(in) :: f
      character(len=*), intent(in) :: scheme
      real(dp), allocatable, intent(out) :: q(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: s(:, :)

      if (trim(scheme) == "chelpg") then
         call chelpg_charges(f%mol, f%density, q, error)
      else
         call f%mol%overlap(s)
         call mulliken_charges(f%mol, f%density, s, q, error)
      end if
   end subroutine fragment_charges

   subroutine report_charges(frag, n_frag, n_atoms, charges, error)
      !! Mulliken charges over the whole system, for the result record
      type(fragment_t), intent(in) :: frag(:)
      integer, intent(in) :: n_frag, n_atoms
      real(dp), allocatable, intent(out) :: charges(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: q(:)
      integer :: f

      allocate (charges(n_atoms), source=0.0_dp)
      do f = 1, n_frag
         call fragment_charges(frag(f), "mulliken", q, error)
         if (error%has_error()) return
         charges(frag(f)%atoms) = q
      end do
   end subroutine report_charges

end module mqc_libcint_fmo
