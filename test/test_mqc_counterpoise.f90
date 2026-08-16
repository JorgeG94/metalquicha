!! Ghost centres, and the basis-set superposition error they exist to expose
!!
!! A counterpoise-corrected many-body term computes each monomer in the *pair's*
!! basis rather than its own. These are the properties that makes rest on, all
!! checkable before any expansion is assembled:
!!
!!   * a ghost centre carries basis functions and no nucleus
!!   * ghosting changes nothing about the AO space -- same count, same ordering
!!   * a monomer in the pair's basis is *lower* than in its own, and the gap is
!!     the BSSE, which is the whole quantity counterpoise removes
!!
!! The last one is the point. Without it a "counterpoise correction" could be
!! wired up, run, and produce a number that is simply the uncorrected one.
module test_mqc_counterpoise
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_counterpoise

   real(dp), parameter :: ANG = 1.8897261254578281_dp

   !> Two waters, 3 Angstrom apart along x. Far enough that the interaction is
   !> small and close enough that each borrows the other's functions.
   real(dp), parameter :: SEP = 3.0_dp

contains

   subroutine collect_counterpoise(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("a_ghost_carries_no_nucleus", test_no_nucleus), &
                  new_unittest("ghosting_preserves_the_ao_space", test_ao_space), &
                  new_unittest("the_pair_basis_lowers_a_monomer", test_bsse_is_real), &
                  new_unittest("an_isolated_ghost_changes_nothing", test_far_ghost) &
                  ]
   end subroutine collect_counterpoise

   subroutine dimer_geometry(z, sym, c)
      !! Two waters, monomer A first then monomer B
      integer, intent(out) :: z(6)
      character(len=2), intent(out) :: sym(6)
      real(dp), intent(out) :: c(3, 6)

      integer :: i

      z = [8, 1, 1, 8, 1, 1]
      sym = ["O ", "H ", "H ", "O ", "H ", "H "]
      c(:, 1) = [0.0_dp, 0.0_dp, 0.10077199_dp]
      c(:, 2) = [0.0_dp, 0.77250895_dp, -0.46780200_dp]
      c(:, 3) = [0.0_dp, -0.77250895_dp, -0.46780200_dp]
      do i = 1, 3
         c(:, i + 3) = c(:, i)
         c(1, i + 3) = c(1, i) + SEP
      end do
      c = c*ANG
   end subroutine dimer_geometry

   subroutine test_no_nucleus(error)
      !! Ghosting monomer B removes its charge and leaves its functions
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: dimer, ghosted
      type(error_t) :: err
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: c(3, 6)
      logical :: ghost(6)

      call dimer_geometry(z, sym, c)
      ghost = [.false., .false., .false., .true., .true., .true.]

      call build_libcint_molecule(z, sym, c, "sto-3g", dimer, err)
      call check(error,.not. err%has_error(), "dimer: "//err%get_full_trace())
      if (allocated(error)) return

      call build_libcint_molecule(z, sym, c, "sto-3g", ghosted, err, ghost=ghost)
      call check(error,.not. err%has_error(), "ghosted: "//err%get_full_trace())
      if (allocated(error)) return

      ! Ten electrons' worth of nucleus gone -- one water -- and nothing else.
      call check(error, nint(sum(dimer%charges)), 20, &
                 "the dimer should carry twenty protons")
      if (allocated(error)) return
      call check(error, nint(sum(ghosted%charges)), 10, &
                 "ghosting one water should leave ten")

      call dimer%destroy()
      call ghosted%destroy()
   end subroutine test_no_nucleus

   subroutine test_ao_space(error)
      !! Same number of basis functions, in the same order
      !!
      !! This is the invariant the whole correction rests on. If ghosting moved
      !! or dropped a function, a monomer-in-pair-basis energy would not be
      !! comparable with the pair's and the difference would be meaningless.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: dimer, ghost_a, ghost_b
      type(error_t) :: err
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: c(3, 6)
      real(dp), allocatable :: s_dimer(:, :), s_ghost(:, :)

      call dimer_geometry(z, sym, c)
      call build_libcint_molecule(z, sym, c, "sto-3g", dimer, err)
      call build_libcint_molecule(z, sym, c, "sto-3g", ghost_a, err, &
                                  ghost=[.true., .true., .true., .false., .false., .false.])
      call build_libcint_molecule(z, sym, c, "sto-3g", ghost_b, err, &
                                  ghost=[.false., .false., .false., .true., .true., .true.])
      call check(error,.not. err%has_error(), "build: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, ghost_a%nao, dimer%nao, &
                 "ghosting A changed the number of basis functions")
      if (allocated(error)) return
      call check(error, ghost_b%nao, dimer%nao, &
                 "ghosting B changed the number of basis functions")
      if (allocated(error)) return

      ! Ordering too, not just the count: the overlap is built from the same
      ! functions in the same places, so it must come back identical.
      call dimer%overlap(s_dimer)
      call ghost_b%overlap(s_ghost)
      call check(error, maxval(abs(s_dimer - s_ghost)), 0.0_dp, &
                 "ghosting changed the overlap matrix, so the AO ordering moved", &
                 thr=1.0e-14_dp)

      call dimer%destroy()
      call ghost_a%destroy()
      call ghost_b%destroy()
   end subroutine test_ao_space

   subroutine test_bsse_is_real(error)
      !! A monomer in the pair's basis lies below the same monomer alone
      !!
      !! The gap is the basis-set superposition error: monomer A, having
      !! borrowed B's functions, describes itself better than its own basis
      !! allows. In a plain expansion that borrowing lands in the pair term and
      !! nowhere else, which is why the dimer looks more bound than it is.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: alone, in_pair
      type(rhf_result_t) :: scf_alone, scf_pair
      type(error_t) :: err
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: c(3, 6), bsse

      call dimer_geometry(z, sym, c)

      ! Monomer A in its own basis: three atoms, nothing else.
      call build_libcint_molecule(z(1:3), sym(1:3), c(:, 1:3), "sto-3g", alone, err)
      call check(error,.not. err%has_error(), "alone: "//err%get_full_trace())
      if (allocated(error)) return
      call run_libcint_rhf(alone, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf_alone, err)
      call check(error,.not. err%has_error(), "alone SCF: "//err%get_full_trace())
      if (allocated(error)) return

      ! The same monomer, same ten electrons, in the pair's basis.
      call build_libcint_molecule(z, sym, c, "sto-3g", in_pair, err, &
                                  ghost=[.false., .false., .false., .true., .true., .true.])
      call check(error,.not. err%has_error(), "in pair: "//err%get_full_trace())
      if (allocated(error)) return
      call run_libcint_rhf(in_pair, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf_pair, err)
      call check(error,.not. err%has_error(), "pair SCF: "//err%get_full_trace())
      if (allocated(error)) return

      bsse = scf_alone%energy - scf_pair%energy

      ! Strictly lower. A variational method handed more functions cannot do
      ! worse, and at 3 Angstrom in a minimal basis it does measurably better.
      call check(error, bsse > 0.0_dp, &
                 "the pair basis did not lower the monomer, so the ghost "// &
                 "functions are not reaching the SCF")
      if (allocated(error)) return

      ! And it is a real effect rather than convergence noise: sto-3g at this
      ! separation is worth more than a microhartree and less than a hartree.
      call check(error, bsse > 1.0e-6_dp, &
                 "the lowering is too small to be superposition error")
      if (allocated(error)) return
      call check(error, bsse < 1.0_dp, &
                 "the lowering is far too large to be superposition error")

      call alone%destroy()
      call in_pair%destroy()
   end subroutine test_bsse_is_real

   subroutine test_far_ghost(error)
      !! Ghost functions a long way off change nothing
      !!
      !! The counterpart to the test above. Superposition error comes from
      !! functions near enough to be borrowed, so at two hundred Angstrom the
      !! ghosted monomer must return to its isolated energy. Without this, a
      !! ghost that was silently ignored and a ghost that worked would look the
      !! same from the sign of one difference.
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: alone, far
      type(rhf_result_t) :: scf_alone, scf_far
      type(error_t) :: err
      integer :: z(6)
      character(len=2) :: sym(6)
      real(dp) :: c(3, 6)
      integer :: i

      call dimer_geometry(z, sym, c)
      call build_libcint_molecule(z(1:3), sym(1:3), c(:, 1:3), "sto-3g", alone, err)
      call run_libcint_rhf(alone, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf_alone, err)
      call check(error,.not. err%has_error(), "alone: "//err%get_full_trace())
      if (allocated(error)) return

      do i = 4, 6
         c(1, i) = c(1, i) + 200.0_dp*ANG
      end do
      call build_libcint_molecule(z, sym, c, "sto-3g", far, err, &
                                  ghost=[.false., .false., .false., .true., .true., .true.])
      call run_libcint_rhf(far, 10, 200, 1.0e-11_dp, 1.0e-9_dp, .false., scf_far, err)
      call check(error,.not. err%has_error(), "far: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, scf_far%energy, scf_alone%energy, &
                 "distant ghost functions changed the energy, which they cannot do", &
                 thr=1.0e-8_dp)

      call alone%destroy()
      call far%destroy()
   end subroutine test_far_ghost

end module test_mqc_counterpoise

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_counterpoise, only: collect_counterpoise
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_counterpoise", collect_counterpoise)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
