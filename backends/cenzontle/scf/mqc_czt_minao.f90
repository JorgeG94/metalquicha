!! The minimal-ANO starting density, PySCF's `minao` guess
module mqc_czt_minao
   !! A starting density built in a minimal atomic natural orbital basis and
   !! projected into the calculation's own.
   !!
   !! Each atom contributes the leading ANO-RCC contractions for its occupied
   !! shells, occupied with the spherically averaged ground-state configuration:
   !! for each l, floor(n_l / 2(2l+1)) doubly occupied columns and then one open
   !! column holding the remainder spread evenly over its 2l+1 components. That
   !! diagonal density is carried into the target basis by `project_density`.
   !!
   !! There is no free-atom SCF, so nothing to converge or cache, and the atoms
   !! are neutral: the guess carries sum(Z) electrons whatever the molecule's
   !! charge, and the first diagonalisation places the rest. PySCF does not
   !! renormalise the electron count either.
   !!
   !! The basis is `basis_sets/minao/minao.json`, H-Kr, written by
   !! `tools/minao/gen_minao_basis.py` from PySCF's `ano` data. The
   !! configurations are PySCF's `NRSRHF_CONFIGURATION`, carried here verbatim.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_elements, only: element_number_to_symbol
   use mqc_string_utils, only: int_to_text
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, atom_ao_blocks, &
                                subshell_layout
   use mqc_czt_projection, only: project_density
   implicit none
   private

   public :: build_minao_guess    !! Target-basis density from the minimal ANO atoms
   public :: minao_occupation     !! Doubly occupied columns and open-shell fraction

   character(len=*), parameter :: MINAO_BASIS = "minao/minao"
      !! Resolved through `find_basis_file`, so `<dir>/minao/minao.json` on the
      !! basis search path. A subdirectory because `basis_sets/*.json` is
      !! derived from `MQC_BASIS_SETS` and anything else there is deleted.

   integer, parameter :: MINAO_Z_MAX = 36

   ! Electrons per angular momentum (s, p, d, f) in the spherically averaged
   ! ground state, from `pyscf.data.elements.NRSRHF_CONFIGURATION` (PySCF
   ! 2.14), copied rather than re-derived so the guess is PySCF's. Its oddities
   ! are kept on purpose: Sc puts its third valence electron in p, Cr is
   ! 4s2 3d4, and Mn-Ni empty 4s.
   integer, parameter :: CONFIGURATION(0:3, MINAO_Z_MAX) = reshape([ &
                                                                   1, 0, 0, 0, &   ! H
                                                                   2, 0, 0, 0, &   ! He
                                                                   3, 0, 0, 0, &   ! Li
                                                                   4, 0, 0, 0, &   ! Be
                                                                   4, 1, 0, 0, &   ! B
                                                                   4, 2, 0, 0, &   ! C
                                                                   4, 3, 0, 0, &   ! N
                                                                   4, 4, 0, 0, &   ! O
                                                                   4, 5, 0, 0, &   ! F
                                                                   4, 6, 0, 0, &   ! Ne
                                                                   5, 6, 0, 0, &   ! Na
                                                                   6, 6, 0, 0, &   ! Mg
                                                                   6, 7, 0, 0, &   ! Al
                                                                   6, 8, 0, 0, &   ! Si
                                                                   6, 9, 0, 0, &   ! P
                                                                   6, 10, 0, 0, &  ! S
                                                                   6, 11, 0, 0, &  ! Cl
                                                                   6, 12, 0, 0, &  ! Ar
                                                                   7, 12, 0, 0, &  ! K
                                                                   8, 12, 0, 0, &  ! Ca
                                                                   8, 13, 0, 0, &  ! Sc
                                                                   8, 12, 2, 0, &  ! Ti
                                                                   8, 12, 3, 0, &  ! V
                                                                   8, 12, 4, 0, &  ! Cr
                                                                   6, 12, 7, 0, &  ! Mn
                                                                   6, 12, 8, 0, &  ! Fe
                                                                   6, 12, 9, 0, &  ! Co
                                                                   6, 12, 10, 0, &  ! Ni
                                                                   7, 12, 10, 0, &  ! Cu
                                                                   8, 12, 10, 0, &  ! Zn
                                                                   8, 13, 10, 0, &  ! Ga
                                                                   8, 14, 10, 0, &  ! Ge
                                                                   8, 15, 10, 0, &  ! As
                                                                   8, 16, 10, 0, &  ! Se
                                                                   8, 17, 10, 0, &  ! Br
                                                                   8, 18, 10, 0 &  ! Kr
                                                                   ], [4, MINAO_Z_MAX])

contains

   pure subroutine minao_occupation(atomic_number, l, n_double, open_fraction)
      !! The occupation pattern of one angular momentum, as PySCF's `frac_occ`
      !!
      !! `n_double` columns hold two electrons per function; the next holds
      !! `open_fraction` per function, which may be zero. Zero and zero for an
      !! element outside H-Kr or an l above f.
      integer, intent(in) :: atomic_number
      integer, intent(in) :: l
      integer, intent(out) :: n_double
      real(dp), intent(out) :: open_fraction   !! Electrons per function, in [0, 2)

      integer :: n_elec, per_shell

      n_double = 0
      open_fraction = 0.0_dp
      if (atomic_number < 1 .or. atomic_number > MINAO_Z_MAX) return
      if (l < 0 .or. l > 3) return
      n_elec = CONFIGURATION(l, atomic_number)
      per_shell = 2*(2*l + 1)
      n_double = n_elec/per_shell
      open_fraction = 2.0_dp*real(n_elec - n_double*per_shell, dp)/real(per_shell, dp)
   end subroutine minao_occupation

   subroutine build_minao_guess(mol, density, error, threshold)
      !! The minao starting density in `mol`'s basis, both spins together
      !!
      !! Refused for an element above Kr, for an atom carrying an effective core
      !! potential, and for a Cartesian basis on an atom with occupied d, whose
      !! six Cartesian functions are not the five the occupations describe.
      !! Ghost centres are given no electrons.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), allocatable, intent(out) :: density(:, :)   !! (nao, nao), total
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: threshold
         !! The SCF's linear-dependence threshold, so the projection drops what
         !! the SCF's orthogonaliser drops. See `project_density`.

      type(czt_molecule_t) :: minimal
      character(len=2), allocatable :: symbols(:)
      logical, allocatable :: ghost(:)
      real(dp), allocatable :: d_small(:, :), occupation(:)
      integer :: iatom, z, n_double, i
      real(dp) :: open_fraction

      if (error%has_error()) return
      if (.not. allocated(mol%atomic_numbers)) then
         call error%set(ERROR_VALIDATION, "minao guess: the molecule does not record "// &
                        "its elements")
         return
      end if
      if (allocated(mol%core_electrons)) then
         if (any(mol%core_electrons > 0)) then
            call error%set(ERROR_VALIDATION, "minao guess: not available with an "// &
                           "effective core potential")
            return
         end if
      end if

      allocate (symbols(mol%natm), ghost(mol%natm))
      do iatom = 1, mol%natm
         z = mol%atomic_numbers(iatom)
         if (z < 1 .or. z > MINAO_Z_MAX) then
            call error%set(ERROR_VALIDATION, "minao guess: the minimal ANO data covers "// &
                           "H-Kr and atom "//int_to_text(iatom)//" has Z="//int_to_text(z))
            return
         end if
         if (mol%cartesian) then
            call minao_occupation(z, 2, n_double, open_fraction)
            if (n_double > 0 .or. open_fraction > 0.0_dp) then
               call error%set(ERROR_VALIDATION, "minao guess: a Cartesian basis on an "// &
                              "element with occupied d (Z="//int_to_text(z)//") is not "// &
                              "supported")
               return
            end if
         end if
         symbols(iatom) = element_number_to_symbol(z)
         ghost(iatom) = mol%charges(iatom) <= 0.0_dp
      end do

      call build_czt_molecule(mol%atomic_numbers, symbols, mol%coords, MINAO_BASIS, &
                              minimal, error, ghost=ghost, force_cartesian=mol%cartesian)
      if (error%has_error()) then
         call error%add_context("minao guess")
         return
      end if

      call minimal_occupations(minimal, ghost, occupation, error)
      if (error%has_error()) then
         call minimal%destroy()
         return
      end if

      allocate (d_small(minimal%nao, minimal%nao))
      d_small = 0.0_dp
      do i = 1, minimal%nao
         d_small(i, i) = occupation(i)
      end do

      call project_density(mol, minimal, d_small, density, error, threshold)
      call minimal%destroy()
   end subroutine build_minao_guess

   subroutine minimal_occupations(minimal, ghost, occupation, error)
      !! Per-function occupation of the minimal molecule
      !!
      !! Walks the contraction columns atom by atom, numbering them within each
      !! l, and checks the file holds exactly the columns the configuration
      !! occupies -- the generator and this table are two copies of one rule.
      type(czt_molecule_t), intent(in) :: minimal
      logical, intent(in) :: ghost(:)
      real(dp), allocatable, intent(out) :: occupation(:)   !! (nao)
      type(error_t), intent(inout) :: error

      integer, allocatable :: ang(:), first(:), ncomp(:), offsets(:), counts(:)
      integer :: seen(0:3)
      integer :: n_sub, a, iatom, z, l, n_double, expected
      real(dp) :: open_fraction, occ

      call subshell_layout(minimal, ang, first, ncomp, n_sub)
      allocate (offsets(minimal%natm), counts(minimal%natm))
      call atom_ao_blocks(minimal, offsets, counts)
      allocate (occupation(minimal%nao))
      occupation = 0.0_dp

      a = 1
      do iatom = 1, minimal%natm
         z = minimal%atomic_numbers(iatom)
         seen = 0
         ! Subshells come in AO order and each atom's functions are one
         ! contiguous run, so this atom's columns are the next ones up to the
         ! end of its block.
         do while (a <= n_sub)
            if (first(a) >= offsets(iatom) + counts(iatom)) exit
            l = ang(a)
            if (l > 3) then
               call error%set(ERROR_VALIDATION, "minao guess: the minimal basis has an "// &
                              "l > 3 shell")
               return
            end if
            seen(l) = seen(l) + 1
            call minao_occupation(z, l, n_double, open_fraction)
            if (seen(l) <= n_double) then
               occ = 2.0_dp
            else
               occ = open_fraction
            end if
            if (.not. ghost(iatom)) occupation(first(a) + 1:first(a) + ncomp(a)) = occ
            a = a + 1
         end do

         do l = 0, 3
            call minao_occupation(z, l, n_double, open_fraction)
            expected = n_double
            if (open_fraction > 0.0_dp) expected = expected + 1
            if (seen(l) /= expected) then
               call error%set(ERROR_VALIDATION, "minao guess: the minimal basis has "// &
                              int_to_text(seen(l))//" columns of l="//int_to_text(l)// &
                              " for Z="//int_to_text(z)//" where the configuration "// &
                              "occupies "//int_to_text(expected))
               return
            end if
         end do
      end do
   end subroutine minimal_occupations

end module mqc_czt_minao
