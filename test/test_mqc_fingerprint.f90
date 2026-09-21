!! Tests for the calculation fingerprint
module test_mqc_fingerprint
   !! The fingerprint has exactly two ways to be wrong, and they are not
   !! equally bad.
   !!
   !! A **false mismatch** -- two identical calculations hashing differently --
   !! costs a recomputation. Annoying, visible, harmless.
   !!
   !! A **false match** -- two different calculations hashing the same -- is
   !! the whole reason this module exists. It means a restart splices energies
   !! from one calculation into another, finishes early, and reports a total
   !! that no check downstream can question. So the bulk of what follows is one
   !! test per material field, each changing that field alone and demanding the
   !! hash move. A field added to a config and not to `add_method` fails here
   !! rather than in someone's results.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp, int32
   use mqc_fingerprint, only: calculation_fingerprint
   use mqc_physical_fragment, only: system_geometry_t, bond_t
   use mqc_method_config, only: method_config_t
   use mqc_method_types, only: METHOD_TYPE_GFN2, METHOD_TYPE_DFT, METHOD_TYPE_HF
   use mqc_calc_types, only: CALC_TYPE_ENERGY, CALC_TYPE_GRADIENT
   implicit none
   private

   public :: collect_mqc_fingerprint

contains

   subroutine collect_mqc_fingerprint(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("deterministic", test_deterministic), &
                  new_unittest("geometry_moves_it", test_geometry), &
                  new_unittest("negative_zero_does_not", test_negative_zero), &
                  new_unittest("elements_move_it", test_elements), &
                  new_unittest("charge_moves_it", test_charge), &
                  new_unittest("partition_moves_it", test_partition), &
                  new_unittest("bonds_move_it", test_bonds), &
                  new_unittest("method_moves_it", test_method), &
                  new_unittest("basis_moves_it", test_basis), &
                  new_unittest("functional_moves_it", test_functional), &
                  new_unittest("scf_threshold_moves_it", test_threshold), &
                  new_unittest("scf_gradient_threshold_moves_it", test_gradient_threshold), &
                  new_unittest("driver_moves_it", test_driver), &
                  new_unittest("verbosity_does_not", test_verbosity), &
                  new_unittest("angular_form_moves_it", test_cartesian), &
                  new_unittest("root_count_moves_it", test_excited_roots), &
                  new_unittest("an_absent_excited_block_does_not", test_excited_absent), &
                  new_unittest("a_fixed_calculation_keeps_its_hash", test_golden_hash) &
                  ]
   end subroutine collect_mqc_fingerprint

   ! -- the two halves of the contract ---------------------------------------

   subroutine test_deterministic(error)
      !! The same calculation, built twice, hashes the same
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: a, b
      type(method_config_t) :: ca, cb

      call water_dimer(a)
      call water_dimer(b)
      call gfn2(ca)
      call gfn2(cb)
      call check(error, calculation_fingerprint(a, ca, CALC_TYPE_ENERGY) == &
                 calculation_fingerprint(b, cb, CALC_TYPE_ENERGY), &
                 "identical calculations must fingerprint identically")
   end subroutine test_deterministic

   subroutine test_verbosity(error)
      !! Things that do not change the number must not change the hash
      !!
      !! A checkpoint written by a quiet batch run has to be usable by a noisy
      !! interactive one. If this fails, the fingerprint is unusable rather
      !! than merely strict.
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(method_config_t) :: quiet, loud

      call water_dimer(sys)
      call gfn2(quiet)
      call gfn2(loud)
      loud%verbose = .true.
      loud%device_rank = 3
      call check(error, calculation_fingerprint(sys, quiet, CALC_TYPE_ENERGY) == &
                 calculation_fingerprint(sys, loud, CALC_TYPE_ENERGY), &
                 "verbosity and GPU binding must not change the fingerprint")
   end subroutine test_verbosity

   subroutine test_negative_zero(error)
      !! An atom at -0.0 is at 0.0
      !!
      !! Same molecule, and the two spellings arrive depending on how the
      !! geometry was built. Without normalising, a checkpoint would be
      !! refused for a difference that does not exist.
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: a, b
      type(method_config_t) :: config

      call water_dimer(a)
      call water_dimer(b)
      b%coordinates(1, 1) = -0.0_dp
      a%coordinates(1, 1) = 0.0_dp
      call gfn2(config)
      call check(error, calculation_fingerprint(a, config, CALC_TYPE_ENERGY) == &
                 calculation_fingerprint(b, config, CALC_TYPE_ENERGY), &
                 "-0.0 and 0.0 are the same position")
   end subroutine test_negative_zero

   ! -- one per material field -----------------------------------------------

   subroutine test_geometry(error)
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: a, b
      type(method_config_t) :: config

      call water_dimer(a)
      call water_dimer(b)
      ! Far below any threshold anyone would notice by eye, and it still has
      ! to move the hash: a resumed run must not splice energies from a
      ! geometry that drifted.
      b%coordinates(1, 1) = b%coordinates(1, 1) + 1.0e-12_dp
      call gfn2(config)
      call differ(error, a, b, config, config, "a 1e-12 Bohr shift")
   end subroutine test_geometry

   subroutine test_elements(error)
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: a, b
      type(method_config_t) :: config

      call water_dimer(a)
      call water_dimer(b)
      b%element_numbers(4) = 9   ! O -> F, same shape, different molecule
      call gfn2(config)
      call differ(error, a, b, config, config, "swapping an element")
   end subroutine test_elements

   subroutine test_charge(error)
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: a, b
      type(method_config_t) :: config

      call water_dimer(a)
      call water_dimer(b)
      b%charge = 1
      call gfn2(config)
      call differ(error, a, b, config, config, "changing the charge")
   end subroutine test_charge

   subroutine test_partition(error)
      !! Same atoms, different grouping
      !!
      !! The one a naive fingerprint misses. Every array has the same shape and
      !! the same contents; only the assignment of atoms to monomers moved --
      !! and with it the meaning of every term index in a stored list.
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: a, b
      type(method_config_t) :: config

      call water_dimer(a)
      call water_dimer(b)
      b%fragment_atoms(1, 1) = 3
      b%fragment_atoms(1, 2) = 0
      call gfn2(config)
      call differ(error, a, b, config, config, "repartitioning the monomers")
   end subroutine test_partition

   subroutine test_bonds(error)
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: a, b
      type(method_config_t) :: config

      call water_dimer(a)
      call water_dimer(b)
      allocate (b%bonds(1))
      b%bonds(1)%atom_i = 0
      b%bonds(1)%atom_j = 3
      b%bonds(1)%order = 1
      b%bonds(1)%is_broken = .false.
      call gfn2(config)
      call differ(error, a, b, config, config, "declaring a bond")
   end subroutine test_bonds

   subroutine test_method(error)
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(method_config_t) :: a, b

      call water_dimer(sys)
      call gfn2(a)
      call gfn2(b)
      b%method_type = METHOD_TYPE_HF
      call differ(error, sys, sys, a, b, "changing the method")
   end subroutine test_method

   subroutine test_basis(error)
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(method_config_t) :: a, b

      call water_dimer(sys)
      call dft(a)
      call dft(b)
      b%basis_set = "def2-tzvp"
      call differ(error, sys, sys, a, b, "changing the basis")
   end subroutine test_basis

   subroutine test_cartesian(error)
      !! `model.cartesian` is a different model, not a different spelling
      !!
      !! Above p the Cartesian and spherical forms span different spaces, so
      !! water at cc-pVDZ is 24 functions and -76.420342 one way and 25 and
      !! -76.421536 the other. Nothing else in the deck distinguishes them --
      !! same basis name, same functional, same thresholds -- so without this
      !! in the hash two genuinely different energies share a fingerprint,
      !! which is the shape a caching bug takes.
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(method_config_t) :: a, b

      call water_dimer(sys)
      call dft(a)
      call dft(b)
      b%scf%cartesian = .true.
      call differ(error, sys, sys, a, b, "forcing a Cartesian basis")
   end subroutine test_cartesian

   subroutine test_functional(error)
      !! The one with no other symptom
      !!
      !! Every dimension, every count, every threshold identical -- only the
      !! functional differs, and the energies are completely different numbers.
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(method_config_t) :: a, b

      call water_dimer(sys)
      call dft(a)
      call dft(b)
      b%dft%functional = "pbe0"
      call differ(error, sys, sys, a, b, "changing the functional")
   end subroutine test_functional

   subroutine test_threshold(error)
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(method_config_t) :: a, b

      call water_dimer(sys)
      call dft(a)
      call dft(b)
      b%scf%energy_convergence = 1.0e-6_dp
      call differ(error, sys, sys, a, b, "loosening the SCF threshold")
   end subroutine test_threshold

   subroutine test_gradient_threshold(error)
      !! `keywords.scf.gradient_tolerance` is half of the convergence test
      !!
      !! Both backends stop on `|dE|` *and* the commutator norm, so a run that
      !! names this one stops at a different iteration and holds a different
      !! density from one that leaves it at the derived default -- with the
      !! same energy threshold, the same basis and the same functional. Nothing
      !! else in the deck distinguishes the two.
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(method_config_t) :: a, b

      call water_dimer(sys)
      call dft(a)
      call dft(b)
      b%scf%gradient_convergence = 1.0e-7_dp
      call differ(error, sys, sys, a, b, "naming a gradient threshold")
   end subroutine test_gradient_threshold

   subroutine test_driver(error)
      !! An energy run and a gradient run are not interchangeable
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(method_config_t) :: config

      call water_dimer(sys)
      call gfn2(config)
      call check(error, calculation_fingerprint(sys, config, CALC_TYPE_ENERGY) /= &
                 calculation_fingerprint(sys, config, CALC_TYPE_GRADIENT), &
                 "energy and gradient runs must not share a fingerprint")
   end subroutine test_driver

   subroutine test_excited_roots(error)
      !! Three roots and ten roots are different calculations
      !!
      !! They share a reference, a basis, a functional and a total energy, so
      !! every other field in the hash agrees. Without the excited block a
      !! restart would splice a three-root checkpoint into a ten-root deck,
      !! finish early, and report seven states it never computed.
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(method_config_t) :: a, b

      call water_dimer(sys)
      call dft(a)
      call dft(b)
      a%excited%enabled = .true.
      a%excited%n_states = 3
      b%excited%enabled = .true.
      b%excited%n_states = 10
      call differ(error, sys, sys, a, b, "asking for more roots")
      if (allocated(error)) return

      ! And the approximation, which changes every energy in the list while
      ! leaving their count alone.
      call dft(a)
      call dft(b)
      a%excited%enabled = .true.
      a%excited%n_states = 5
      a%excited%method = "tda"
      b%excited%enabled = .true.
      b%excited%n_states = 5
      b%excited%method = "rpa"
      call differ(error, sys, sys, a, b, "switching TDA for RPA")
      if (allocated(error)) return

      call dft(a)
      call dft(b)
      a%excited%enabled = .true.
      a%excited%n_states = 5
      a%excited%spin = "singlet"
      b%excited%enabled = .true.
      b%excited%n_states = 5
      b%excited%spin = "triplet"
      call differ(error, sys, sys, a, b, "asking for triplets instead of singlets")
   end subroutine test_excited_roots

   subroutine test_excited_absent(error)
      !! No excited states asked for hashes as it always did
      !!
      !! The other half of the Layer 1 contract: a deck with no
      !! `keywords.excited_states` must produce the same fingerprint it
      !! produced before the block existed, or every checkpoint ever written
      !! is invalidated by a feature nobody used. The block's own defaults are
      !! deliberately not hashed when it is off, which is what this pins.
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(method_config_t) :: a, b

      call water_dimer(sys)
      call dft(a)
      call dft(b)
      ! `b` carries settings a deck could have left lying around with the
      ! block switched off. None of them may reach the hash.
      b%excited%n_states = 0
      b%excited%method = "tda"
      b%excited%spin = "triplet"
      b%excited%tolerance = 1.0e-7_dp
      call check(error, calculation_fingerprint(sys, a, CALC_TYPE_ENERGY) == &
                 calculation_fingerprint(sys, b, CALC_TYPE_ENERGY), &
                 "an excited block that asks for no roots must not move the hash")
   end subroutine test_excited_absent

   subroutine test_golden_hash(error)
      !! One fingerprint, written down, for a calculation that never changes
      !!
      !! Every other test here is relative -- this field moved, so the hash
      !! must move -- and a relative test is blind to the hash moving under a
      !! calculation nobody touched. That is the expensive failure: a field
      !! hashed unconditionally, outside the `if (config%excited%enabled)`
      !! guard or any other, invalidates every checkpoint ever written and
      !! passes every test above. A literal is the only thing that catches it.
      !!
      !! When this fails and the change was deliberate, move
      !! `mqc-fingerprint-v1` in `calculation_fingerprint` with it, so stored
      !! energies are rejected wholesale rather than mismatched one at a time,
      !! and then write the new digest down here.
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t) :: sys
      type(method_config_t) :: config

      character(len=*), parameter :: GOLDEN = "0684a824686c71d4"

      call water_dimer(sys)
      call dft(config)
      call check(error, calculation_fingerprint(sys, config, CALC_TYPE_ENERGY) == GOLDEN, &
                 "b3lyp/def2-svp on the water dimer used to fingerprint as "//GOLDEN// &
                 " and now hashes as "//calculation_fingerprint(sys, config, CALC_TYPE_ENERGY)// &
                 ". Every checkpoint written before this change is now unreadable")
   end subroutine test_golden_hash

   ! -- helpers ---------------------------------------------------------------

   subroutine differ(error, a, b, ca, cb, what)
      !! Demand that a change moves the fingerprint
      type(error_type), allocatable, intent(out) :: error
      type(system_geometry_t), intent(in) :: a, b
      type(method_config_t), intent(in) :: ca, cb
      character(len=*), intent(in) :: what

      call check(error, calculation_fingerprint(a, ca, CALC_TYPE_ENERGY) /= &
                 calculation_fingerprint(b, cb, CALC_TYPE_ENERGY), &
                 what//" must change the fingerprint, or a restart will reuse "// &
                 "energies across it")
   end subroutine differ

   subroutine water_dimer(sys)
      !! Two waters, six atoms, two monomers. Bohr, 0-based atom indices.
      type(system_geometry_t), intent(out) :: sys

      integer :: i

      sys%n_monomers = 2
      sys%atoms_per_monomer = 3
      sys%total_atoms = 6
      sys%charge = 0
      sys%multiplicity = 1

      allocate (sys%element_numbers(6))
      sys%element_numbers = [8, 1, 1, 8, 1, 1]
      allocate (sys%coordinates(3, 6))
      do i = 1, 6
         sys%coordinates(:, i) = [real(i, dp), 0.5_dp*i, -0.25_dp*i]
      end do

      allocate (sys%fragment_sizes(2))
      sys%fragment_sizes = 3
      allocate (sys%fragment_atoms(3, 2))
      sys%fragment_atoms(:, 1) = [0, 1, 2]
      sys%fragment_atoms(:, 2) = [3, 4, 5]
      allocate (sys%fragment_charges(2))
      sys%fragment_charges = 0
      allocate (sys%fragment_multiplicities(2))
      sys%fragment_multiplicities = 1
   end subroutine water_dimer

   subroutine gfn2(config)
      type(method_config_t), intent(out) :: config
      config%method_type = METHOD_TYPE_GFN2
   end subroutine gfn2

   subroutine dft(config)
      type(method_config_t), intent(out) :: config
      config%method_type = METHOD_TYPE_DFT
      config%basis_set = "def2-svp"
      config%dft%functional = "b3lyp"
   end subroutine dft

end module test_mqc_fingerprint

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_fingerprint, only: collect_mqc_fingerprint
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_fingerprint", collect_mqc_fingerprint)]
   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
