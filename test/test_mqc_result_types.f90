module test_mqc_result_types
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_result_types, only: mp2_energy_t, cc_energy_t, energy_t, calculation_result_t, &
                               mbe_result_t, STATE_SPIN_SINGLET
   use pic_types, only: dp
   use mqc_quao_rows, only: quao_rows_t, QUAO_ROW_BOND, QUAO_ROW_DELOCALIZATION
   use pic_test_helpers, only: is_equal
   implicit none
   private
   public :: collect_mqc_result_types_tests

contains

   !! Collect all exported unit tests
   subroutine collect_mqc_result_types_tests(testsuite)
      !! Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("mp2_energy_total", test_mp2_total), &
                  new_unittest("mp2_energy_scs", test_mp2_scs), &
                  new_unittest("mp2_energy_reset", test_mp2_reset), &
                  new_unittest("cc_energy_total", test_cc_total), &
                  new_unittest("cc_energy_reset", test_cc_reset), &
                  new_unittest("energy_total", test_energy_total), &
                  new_unittest("energy_reset", test_energy_reset), &
                  new_unittest("result_initialization", test_result_init), &
                  new_unittest("result_destroy", test_result_destroy), &
                  new_unittest("result_reset", test_result_reset), &
                  new_unittest("quao_rows_pack_round_trip", test_quao_rows_round_trip), &
                  new_unittest("mp2_stability_warns_on_positive", test_mp2_stability), &
                  new_unittest("cc_stability_warns_on_positive", test_cc_stability), &
                  new_unittest("mbe_result_allocates_derivatives", test_mbe_allocate) &
                  ]
   end subroutine collect_mqc_result_types_tests

   subroutine test_mp2_total(error)
      type(error_type), allocatable, intent(out) :: error
      type(mp2_energy_t) :: mp2
      real(dp) :: total

      ! Set MP2 components (negative = correlation energy)
      mp2%ss = -0.1_dp
      mp2%os = -0.3_dp

      ! Calculate total
      total = mp2%total()

      ! Check total = ss + os
      call check(error, is_equal(total, -0.4_dp), &
                 "MP2 total should be sum of ss and os")
   end subroutine test_mp2_total

   subroutine test_mp2_scs(error)
      type(error_type), allocatable, intent(out) :: error
      type(mp2_energy_t) :: mp2
      real(dp) :: scs_energy, expected_scs

      ! Set MP2 components (negative = correlation energy)
      mp2%ss = -0.3_dp
      mp2%os = -0.6_dp

      ! Calculate SCS-MP2: (1/3)*ss + 1.2*os
      scs_energy = mp2%scs()
      expected_scs = (1.0_dp/3.0_dp)*(-0.3_dp) + 1.2_dp*(-0.6_dp)

      ! Check SCS calculation
      call check(error, is_equal(scs_energy, expected_scs), &
                 "SCS-MP2 should use correct scaling factors")
      if (allocated(error)) return

      ! Verify expected value explicitly: -0.1 + (-0.72) = -0.82
      call check(error, is_equal(scs_energy, -0.82_dp), &
                 "SCS-MP2 value should be -0.82")
   end subroutine test_mp2_scs

   subroutine test_mp2_reset(error)
      type(error_type), allocatable, intent(out) :: error
      type(mp2_energy_t) :: mp2

      ! Set non-zero values
      mp2%ss = -1.5_dp
      mp2%os = -2.5_dp

      ! Reset
      call mp2%reset()

      ! Check all components are zero
      call check(error, is_equal(mp2%ss, 0.0_dp), "MP2 ss should be zero after reset")
      if (allocated(error)) return

      call check(error, is_equal(mp2%os, 0.0_dp), "MP2 os should be zero after reset")
   end subroutine test_mp2_reset

   subroutine test_cc_total(error)
      type(error_type), allocatable, intent(out) :: error
      type(cc_energy_t) :: cc
      real(dp) :: total

      ! Set CC components (negative = correlation energy)
      cc%singles = -0.05_dp
      cc%doubles = -0.25_dp
      cc%triples = -0.10_dp

      ! Calculate total
      total = cc%total()

      ! Check total = singles + doubles + triples
      call check(error, is_equal(total, -0.40_dp), &
                 "CC total should be sum of all components")
   end subroutine test_cc_total

   subroutine test_cc_reset(error)
      type(error_type), allocatable, intent(out) :: error
      type(cc_energy_t) :: cc

      ! Set non-zero values
      cc%singles = -1.0_dp
      cc%doubles = -2.0_dp
      cc%triples = -3.0_dp

      ! Reset
      call cc%reset()

      ! Check all components are zero
      call check(error, is_equal(cc%singles, 0.0_dp), "CC singles should be zero after reset")
      if (allocated(error)) return

      call check(error, is_equal(cc%doubles, 0.0_dp), "CC doubles should be zero after reset")
      if (allocated(error)) return

      call check(error, is_equal(cc%triples, 0.0_dp), "CC triples should be zero after reset")
   end subroutine test_cc_reset

   subroutine test_energy_total(error)
      type(error_type), allocatable, intent(out) :: error
      type(energy_t) :: energy
      real(dp) :: total, expected

      ! Set all energy components (correlations are negative)
      energy%scf = -10.0_dp
      energy%mp2%ss = -0.1_dp
      energy%mp2%os = -0.2_dp
      energy%cc%singles = -0.05_dp
      energy%cc%doubles = -0.15_dp
      energy%cc%triples = -0.03_dp

      ! Calculate total
      total = energy%total()

      ! Expected: scf + mp2%total() + cc%total()
      ! = -10.0 + (-0.1 + -0.2) + (-0.05 + -0.15 + -0.03)
      ! = -10.0 - 0.3 - 0.23 = -10.53
      expected = -10.53_dp

      call check(error, is_equal(total, expected), &
                 "Total energy should be sum of all components")
   end subroutine test_energy_total

   subroutine test_energy_reset(error)
      type(error_type), allocatable, intent(out) :: error
      type(energy_t) :: energy

      ! Set non-zero values
      energy%scf = -100.0_dp
      energy%mp2%ss = -1.0_dp
      energy%mp2%os = -2.0_dp
      energy%cc%singles = -0.5_dp
      energy%cc%doubles = -1.5_dp
      energy%cc%triples = -0.3_dp

      ! Reset
      call energy%reset()

      ! Check all components are zero
      call check(error, is_equal(energy%scf, 0.0_dp), "SCF energy should be zero after reset")
      if (allocated(error)) return

      call check(error, is_equal(energy%mp2%ss, 0.0_dp), "MP2 ss should be zero after reset")
      if (allocated(error)) return

      call check(error, is_equal(energy%mp2%os, 0.0_dp), "MP2 os should be zero after reset")
      if (allocated(error)) return

      call check(error, is_equal(energy%cc%singles, 0.0_dp), "CC singles should be zero after reset")
      if (allocated(error)) return

      call check(error, is_equal(energy%cc%doubles, 0.0_dp), "CC doubles should be zero after reset")
      if (allocated(error)) return

      call check(error, is_equal(energy%cc%triples, 0.0_dp), "CC triples should be zero after reset")
   end subroutine test_energy_reset

   subroutine test_result_init(error)
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      ! Check initial state
      call check(error, is_equal(result%energy%total(), 0.0_dp), &
                 "Initial total energy should be zero")
      if (allocated(error)) return

      call check(error,.not. result%has_energy, "has_energy should be false initially")
      if (allocated(error)) return

      call check(error,.not. result%has_gradient, "has_gradient should be false initially")
      if (allocated(error)) return

      call check(error,.not. allocated(result%gradient), "gradient should not be allocated initially")
   end subroutine test_result_init

   subroutine test_result_destroy(error)
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      ! Allocate and populate. Every allocatable the container owns, not just
      ! the gradient: `result_recv` calls destroy precisely so a container can
      ! be reused, so anything left behind comes back on the next fragment at
      ! the previous fragment's size. Bond orders and the three Fukui arrays
      ! were the ones that used to survive.
      allocate (result%gradient(3, 5))
      result%gradient = 1.0_dp
      result%has_gradient = .true.
      result%energy%scf = -50.0_dp

      allocate (result%sigma(3, 3), result%hessian(15, 15))
      result%sigma = 0.01_dp
      result%hessian = 0.02_dp
      allocate (result%dipole(3), result%dipole_derivatives(3, 15))
      result%dipole = 0.03_dp
      result%dipole_derivatives = 0.04_dp

      allocate (result%bond_orders(5, 5))
      result%bond_orders = 0.5_dp
      allocate (result%bond_order_valences(5))
      result%bond_order_valences = 2.0_dp
      allocate (result%atomic_charges(5), result%spin_populations(5))
      result%atomic_charges = 0.05_dp
      result%spin_populations = 0.06_dp
      allocate (result%ieda_atom(5), result%ieda_free_atom(5))
      result%ieda_atom = -1.0_dp
      result%ieda_free_atom = -0.9_dp
      allocate (result%ieda_pair(5, 5), result%ieda_classical(5, 5))
      result%ieda_pair = -0.1_dp
      result%ieda_classical = -0.05_dp
      allocate (result%fukui_plus(5))
      result%fukui_plus = 0.1_dp
      allocate (result%fukui_minus(5))
      result%fukui_minus = 0.2_dp
      allocate (result%fukui_dual(5))
      result%fukui_dual = -0.1_dp

      ! The excited-state block, for the same reason: `result_recv` destroys a
      ! container to reuse it, and pic-mpi's array receives allocate only when
      ! the target is unallocated. A spectrum of three roots surviving into a
      ! fragment with ten would be received into the smaller allocation.
      allocate (result%excitation_energies(3))
      result%excitation_energies = 0.3_dp
      allocate (result%excited_total_energies(3))
      result%excited_total_energies = -49.7_dp
      allocate (result%oscillator_strengths(3))
      result%oscillator_strengths = 0.02_dp
      allocate (result%oscillator_strengths_velocity(3))
      result%oscillator_strengths_velocity = 0.03_dp
      allocate (result%transition_dipoles(3, 3))
      result%transition_dipoles = 0.5_dp
      allocate (result%transition_velocities(3, 3))
      result%transition_velocities = 0.4_dp
      allocate (result%transition_dipole_origin(3))
      result%transition_dipole_origin = 0.7_dp
      allocate (result%nto_leading_weight(3))
      result%nto_leading_weight = 0.95_dp
      allocate (result%state_spin(3))
      result%state_spin = STATE_SPIN_SINGLET
      result%has_excited_states = .true.

      ! The bonding-analysis rows, one of each kind, and their atom-pair sums.
      call result%quao_rows%unpack([QUAO_ROW_BOND, 3, 4, 1, 2, 3, 3, 1, 1, 2, 1], &
                                   [1.0_dp, 1.0_dp, 0.96_dp, -67.8_dp])
      allocate (result%quao_rows%atom_bond_index(5, 5), &
                result%quao_rows%atom_kinetic_bond_order(5, 5))
      result%quao_rows%atom_bond_index = 0.1_dp
      result%quao_rows%atom_kinetic_bond_order = -1.0_dp
      result%has_quao_rows = .true.

      ! Everything the type owns is now allocated. Checked rather than
      ! assumed: a field added to the type and to `mismatch` but not to the
      ! block above would leave this case passing over a field it never
      ! allocated, which is the same blind spot in a different place.
      call check(error, len_trim(mismatch(result, .true.)) == 0, "the destroy case "// &
                 "never allocated "//trim(mismatch(result, .true.)))
      if (allocated(error)) return

      ! Destroy
      call result%destroy()

      ! Check deallocation, field by field and by name.
      call check(error, len_trim(mismatch(result, .false.)) == 0, "result_destroy "// &
                 "left "//trim(mismatch(result, .false.))//" allocated")
      if (allocated(error)) return
      call check(error,.not. result%has_excited_states, &
                 "has_excited_states should be false after destroy")
      if (allocated(error)) return
      call check(error,.not. result%has_quao_rows .and. result%quao_rows%n == 0, &
                 "the bonding rows should be empty after destroy")
      if (allocated(error)) return

      ! Check reset was called
      call check(error,.not. result%has_gradient, &
                 "has_gradient should be false after destroy")
   end subroutine test_result_destroy

   subroutine test_quao_rows_round_trip(error)
      !! `pack_integers`/`pack_reals` and `unpack` are inverses, column by column
      !!
      !! They are the MPI wire format for the bonding rows, so a column packed
      !! into the wrong slot would reach the coordinator as a different atom or
      !! orbital type and nothing downstream could tell.
      type(error_type), allocatable, intent(out) :: error
      type(quao_rows_t) :: sent, received

      call sent%unpack([QUAO_ROW_BOND, 11, 15, 1, 2, 3, 3, 1, 1, 2, 1, &
                        QUAO_ROW_DELOCALIZATION, 13, 9, 2, 1, 1, 3, 1, 1, 0, 3], &
                       [0.99_dp, 1.01_dp, 0.9658_dp, -67.81_dp, &
                        1.93_dp, 1.02_dp, 0.3780_dp, -9.21_dp])
      call check(error, sent%n, 2, "two rows should have been unpacked")
      if (allocated(error)) return
      call check(error, sent%atom(2, 2), 1, "the acceptor atom of row 2")
      if (allocated(error)) return
      call check(error, sent%partner_atom(1, 2), 0, "a lone pair has no partner")
      if (allocated(error)) return

      call received%unpack(sent%pack_integers(), sent%pack_reals())
      call check(error, received%n, sent%n, "row count")
      if (allocated(error)) return
      call check(error, all(received%kind == sent%kind), "kind")
      if (allocated(error)) return
      call check(error, all(received%orbital == sent%orbital), "orbital")
      if (allocated(error)) return
      call check(error, all(received%atom == sent%atom), "atom")
      if (allocated(error)) return
      call check(error, all(received%orbital_type == sent%orbital_type), "orbital_type")
      if (allocated(error)) return
      call check(error, all(received%dominant_l == sent%dominant_l), "dominant_l")
      if (allocated(error)) return
      call check(error, all(received%partner_atom == sent%partner_atom), "partner_atom")
      if (allocated(error)) return
      call check(error, all(received%occupation == sent%occupation), "occupation")
      if (allocated(error)) return
      call check(error, all(received%bond_order == sent%bond_order), "bond_order")
      if (allocated(error)) return
      call check(error, all(received%kinetic_bond_order == sent%kinetic_bond_order), &
                 "kinetic_bond_order")
   end subroutine test_quao_rows_round_trip

   function mismatch(result, want) result(name)
      !! The first component whose allocation is not `want`, by name, or empty
      !!
      !! Every allocatable component `calculation_result_t` declares is
      !! listed, in declaration order, and the list is walked twice: once
      !! with `want = .true.`, which says the case populated everything the
      !! type owns, and once after `destroy` with `want = .false.`, which
      !! says nothing survived. A field added to the type and to
      !! `result_destroy` but not here is the one way past this case, so the
      !! list is the inventory the test is really checking.
      type(calculation_result_t), intent(in) :: result
      logical, intent(in) :: want
      character(len=32) :: name

      name = ""
      call note(allocated(result%gradient), "gradient", want, name)
      call note(allocated(result%sigma), "sigma", want, name)
      call note(allocated(result%hessian), "hessian", want, name)
      call note(allocated(result%dipole), "dipole", want, name)
      call note(allocated(result%dipole_derivatives), "dipole_derivatives", want, name)
      call note(allocated(result%bond_orders), "bond_orders", want, name)
      call note(allocated(result%bond_order_valences), "bond_order_valences", want, name)
      call note(allocated(result%atomic_charges), "atomic_charges", want, name)
      call note(allocated(result%spin_populations), "spin_populations", want, name)
      call note(allocated(result%ieda_atom), "ieda_atom", want, name)
      call note(allocated(result%ieda_free_atom), "ieda_free_atom", want, name)
      call note(allocated(result%ieda_pair), "ieda_pair", want, name)
      call note(allocated(result%ieda_classical), "ieda_classical", want, name)
      call note(allocated(result%fukui_plus), "fukui_plus", want, name)
      call note(allocated(result%fukui_minus), "fukui_minus", want, name)
      call note(allocated(result%fukui_dual), "fukui_dual", want, name)
      call note(allocated(result%excitation_energies), "excitation_energies", want, name)
      call note(allocated(result%excited_total_energies), "excited_total_energies", want, name)
      call note(allocated(result%oscillator_strengths), "oscillator_strengths", want, name)
      call note(allocated(result%oscillator_strengths_velocity), "oscillator_strengths_velocity", want, name)
      call note(allocated(result%transition_dipoles), "transition_dipoles", want, name)
      call note(allocated(result%transition_velocities), "transition_velocities", want, name)
      call note(allocated(result%transition_dipole_origin), "transition_dipole_origin", want, name)
      call note(allocated(result%nto_leading_weight), "nto_leading_weight", want, name)
      call note(allocated(result%state_spin), "state_spin", want, name)
      call note(allocated(result%quao_rows%kind), "quao_rows%kind", want, name)
      call note(allocated(result%quao_rows%orbital), "quao_rows%orbital", want, name)
      call note(allocated(result%quao_rows%atom), "quao_rows%atom", want, name)
      call note(allocated(result%quao_rows%orbital_type), "quao_rows%orbital_type", want, name)
      call note(allocated(result%quao_rows%dominant_l), "quao_rows%dominant_l", want, name)
      call note(allocated(result%quao_rows%partner_atom), "quao_rows%partner_atom", want, name)
      call note(allocated(result%quao_rows%occupation), "quao_rows%occupation", want, name)
      call note(allocated(result%quao_rows%bond_order), "quao_rows%bond_order", want, name)
      call note(allocated(result%quao_rows%kinetic_bond_order), &
                "quao_rows%kinetic_bond_order", want, name)
      call note(allocated(result%quao_rows%atom_bond_index), "quao_rows%atom_bond_index", &
                want, name)
      call note(allocated(result%quao_rows%atom_kinetic_bond_order), &
                "quao_rows%atom_kinetic_bond_order", want, name)
   end function mismatch

   subroutine note(state, field, want, name)
      !! Record `field` as the first mismatch, if it is one and none is yet
      logical, intent(in) :: state, want
      character(len=*), intent(in) :: field
      character(len=*), intent(inout) :: name

      if (len_trim(name) > 0) return
      if (state .eqv. want) return
      name = field
   end subroutine note

   subroutine test_result_reset(error)
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t) :: result

      ! Set various values
      result%energy%scf = -100.0_dp
      result%energy%mp2%ss = -0.5_dp
      result%has_energy = .true.
      result%has_gradient = .true.
      result%has_excited_states = .true.

      ! Reset
      call result%reset()

      ! Check all flags are false
      call check(error,.not. result%has_energy, "has_energy should be false after reset")
      if (allocated(error)) return

      call check(error,.not. result%has_gradient, "has_gradient should be false after reset")
      if (allocated(error)) return

      call check(error,.not. result%has_hessian, "has_hessian should be false after reset")
      if (allocated(error)) return

      call check(error,.not. result%has_dipole, "has_dipole should be false after reset")
      if (allocated(error)) return

      call check(error,.not. result%has_excited_states, &
                 "has_excited_states should be false after reset")
      if (allocated(error)) return

      ! Check energy was reset
      call check(error, is_equal(result%energy%total(), 0.0_dp), &
                 "Total energy should be zero after reset")
   end subroutine test_result_reset

   subroutine test_mp2_stability(error)
      !! A positive MP2 correlation component means something has gone wrong --
      !! the second-order correction is negative-definite for a closed-shell
      !! reference -- so the components are checked and the run warns. What is
      !! asserted here is that asking does not disturb them: a check that
      !! silently corrected the sign would hide the very thing it exists to
      !! report.
      type(error_type), allocatable, intent(out) :: error

      type(mp2_energy_t) :: energy

      energy%ss = 0.01_dp      ! both wrong on purpose
      energy%os = 0.02_dp
      call energy%check_stability()
      call check(error, abs(energy%ss - 0.01_dp) < 1.0e-14_dp, &
                 "the stability check changed the same-spin component")
      if (allocated(error)) return
      call check(error, abs(energy%os - 0.02_dp) < 1.0e-14_dp, &
                 "the stability check changed the opposite-spin component")
      if (allocated(error)) return

      ! And the ordinary case runs clean.
      energy%ss = -0.01_dp
      energy%os = -0.03_dp
      call energy%check_stability()
      call check(error, abs(energy%total() + 0.04_dp) < 1.0e-14_dp, &
                 "a healthy MP2 energy did not survive its own check")
   end subroutine test_mp2_stability

   subroutine test_cc_stability(error)
      !! The same for coupled cluster, where all three components are checked
      !! and the triples are the one that most often comes back positive when
      !! an index is wrong.
      type(error_type), allocatable, intent(out) :: error

      type(cc_energy_t) :: energy

      energy%singles = 0.001_dp
      energy%doubles = 0.002_dp
      energy%triples = 0.003_dp
      call energy%check_stability()
      call check(error, abs(energy%singles - 0.001_dp) < 1.0e-14_dp .and. &
                 abs(energy%doubles - 0.002_dp) < 1.0e-14_dp .and. &
                 abs(energy%triples - 0.003_dp) < 1.0e-14_dp, &
                 "the stability check altered the coupled-cluster components")
      if (allocated(error)) return

      energy%singles = 0.0_dp
      energy%doubles = -0.05_dp
      energy%triples = -0.001_dp
      call energy%check_stability()
      call check(error, abs(energy%total() + 0.051_dp) < 1.0e-14_dp, &
                 "a healthy coupled-cluster energy did not survive its own check")
   end subroutine test_cc_stability

   subroutine test_mbe_allocate(error)
      !! The expansion's own result allocates its derivatives by atom count,
      !! and both are shaped from the same number: a gradient is (3, N) and a
      !! Hessian (3N, 3N). Allocating twice must not leave the first size
      !! behind, which is what the deallocate inside each is for -- a
      !! reallocation onto a different molecule is exactly what a fragmented
      !! run does, once per fragment.
      type(error_type), allocatable, intent(out) :: error

      type(mbe_result_t) :: result

      call result%allocate_gradient(3)
      call check(error, allocated(result%gradient), "no gradient was allocated")
      if (allocated(error)) return
      call check(error, size(result%gradient, 1) == 3 .and. size(result%gradient, 2) == 3, &
                 "the gradient is not (3, natoms)")
      if (allocated(error)) return
      call check(error, all(abs(result%gradient) < 1.0e-14_dp), &
                 "the gradient was not zeroed")
      if (allocated(error)) return

      call result%allocate_hessian(3)
      call check(error, allocated(result%hessian), "no Hessian was allocated")
      if (allocated(error)) return
      call check(error, size(result%hessian, 1) == 9 .and. size(result%hessian, 2) == 9, &
                 "the Hessian is not (3N, 3N)")
      if (allocated(error)) return

      ! Reallocate smaller: the previous shape must be gone, not merely unused.
      call result%allocate_gradient(2)
      call result%allocate_hessian(2)
      call check(error, size(result%gradient, 2) == 2, "the gradient kept its old size")
      if (allocated(error)) return
      call check(error, size(result%hessian, 1) == 6, "the Hessian kept its old size")
      if (allocated(error)) return

      call result%destroy()
      call check(error,.not. allocated(result%gradient) .and. &
                 .not. allocated(result%hessian), "destroy left arrays behind")
   end subroutine test_mbe_allocate

end module test_mqc_result_types

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_result_types, only: collect_mqc_result_types_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_result_types", collect_mqc_result_types_tests) &
                ]

   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if

end program tester
