module test_mqc_czt_neo
   !! NEO-HF: quantum protons against PySCF-NEO
   !!
   !! The reference is the `pyscf/neo` module of Yang Yang's PySCF fork, whose
   !! `test_hf.py` pins HCN with the proton quantised. The two codes share
   !! nothing: PySCF builds the cross Coulomb terms from its own two-basis
   !! integrals and runs one DIIS over all components, this code makes the
   !! cross terms out of a combined molecule and macro-iterates. Agreement to
   !! the macro-iteration's tolerance is therefore agreement on the physics --
   !! the mass, the sign of every coupling, the ghosting of the quantum
   !! nucleus and the block bookkeeping in the combined basis.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_neo, only: neo_result_t, run_czt_neo_hf
   implicit none
   private
   public :: collect_mqc_czt_neo_tests

   real(dp), parameter :: ANG = 1.8897261254578281_dp

contains

   subroutine collect_mqc_czt_neo_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("hcn_with_a_quantum_proton_matches_pyscf_neo", test_hcn), &
                  new_unittest("the_proton_density_is_one_particle", test_normalisation), &
                  new_unittest("a_heavy_atom_is_refused_for_now", test_refusal) &
                  ]
   end subroutine collect_mqc_czt_neo_tests

   subroutine hcn(result, err, quantum)
      !! HCN as PySCF-NEO's test has it: H at the origin, C and N on z, in Angstrom
      type(neo_result_t), intent(out) :: result
      type(error_t), intent(inout) :: err
      logical, intent(in) :: quantum(3)
      real(dp) :: c(3, 3)
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 1.064_dp*ANG, &
                   0.0_dp, 0.0_dp, 2.220_dp*ANG], [3, 3])
      call run_czt_neo_hf([1, 6, 7], ["H ", "C ", "N "], c, "cc-pvdz", "pb4-d", quantum, &
                          14, 200, 1.0e-10_dp, 1.0e-8_dp, .false., result, err, in_core=.true.)
   end subroutine hcn

   subroutine test_hcn(error)
      !! The energy PySCF-NEO gives for the same molecule, basis and proton basis
      !!
      !! Two references, because the electronic basis can be built either way
      !! and the number depends on it: `-92.8437063565785` is the fork's own
      !! test value in spherical harmonics, `-92.8442210525` the same run with
      !! `cart=True`. Both from the fork at commit f9c0266 with PB4-D.
      type(error_type), allocatable, intent(out) :: error
      type(neo_result_t) :: result
      type(error_t) :: err
      real(dp) :: reference

      call hcn(result, err, [.true., .false., .false.])
      call check(error,.not. err%has_error(), "NEO-HF failed: "//err%get_message())
      if (allocated(error)) return
      call check(error, result%converged, "the macro-iteration did not converge")
      if (allocated(error)) return
      if (result%cartesian) then
         reference = -92.8442210525_dp
      else
         reference = -92.8437063565785_dp
      end if
      call check(error, abs(result%energy - reference) < 2.0e-6_dp, &
                 "the NEO-HF energy disagrees with PySCF-NEO")
      if (allocated(error)) return
      ! The proton sits in a well: its orbital energy is negative, and the
      ! electron-proton attraction is what holds it there.
      call check(error, result%nuclear_orbital_energies(1, 1) < 0.0_dp, &
                 "the proton orbital energy is not negative")
      if (allocated(error)) return
      call check(error, result%electron_nucleus < 0.0_dp, &
                 "the electron-proton energy is not attractive")
   end subroutine test_hcn

   subroutine test_normalisation(error)
      !! `Tr(D_p S_p) = 1`: one proton, whatever the basis does
      type(error_type), allocatable, intent(out) :: error
      type(neo_result_t) :: result
      type(error_t) :: err
      real(dp) :: trace

      call hcn(result, err, [.true., .false., .false.])
      call check(error,.not. err%has_error(), "NEO-HF failed: "//err%get_message())
      if (allocated(error)) return
      trace = sum(result%nuclear_densities(:, :, 1)*result%nuclear_overlaps(:, :, 1))
      call check(error, abs(trace - 1.0_dp) < 1.0e-10_dp, "the proton density is not one particle")
   end subroutine test_normalisation

   subroutine test_refusal(error)
      !! Quantising carbon is not implemented yet and must say so, not run
      type(error_type), allocatable, intent(out) :: error
      type(neo_result_t) :: result
      type(error_t) :: err

      call hcn(result, err, [.true., .true., .false.])
      call check(error, err%has_error(), "a quantum carbon was accepted")
   end subroutine test_refusal

end module test_mqc_czt_neo

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_neo, only: collect_mqc_czt_neo_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_neo", collect_mqc_czt_neo_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
