!! That a written fragment potential reads back as the parameters that went in
module test_mqc_efp_read
   !! The writer and the reader are two halves of one format, and the format is
   !! GAMESS's rather than one chosen here -- labels, fixed columns, records
   !! continued with a trailing `>`. That makes a round trip the test worth having:
   !! compute a potential, write it, read it back, and require the parameters to
   !! survive.
   !!
   !! **Why a round trip rather than a committed reference file.** A stored `.efp`
   !! to parse would test the reader against a snapshot, and the snapshot would
   !! quietly stop matching the writer the first time either changed. Going through
   !! both means neither can drift without this failing.
   !!
   !! It does mean the two could agree on something GAMESS would reject, which a
   !! round trip cannot see. That is what the GAMESS comparison in
   !! `tools/efp_validation/dimer_energy.py` is for, and it needs GAMESS installed,
   !! so it cannot live here.
   !!
   !! The tolerance is the writer's own precision, not the arithmetic's: the file
   !! carries ten decimal places, so a value that survives to 1e-9 survives exactly
   !! as far as the format allows.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_efp_potential, only: efp_potential_t, make_efp_potential, &
                                write_efp_potential
   use mqc_efp_read, only: efp_fragment_t, read_efp_potential
   use mqc_error, only: error_t
   use mqc_cgto, only: molecular_basis_type
   use mqc_basis_utils, only: find_basis_file
   use mqc_json_basis_reader, only: build_molecular_basis_json
   implicit none
   private

   public :: collect_mqc_efp_read_tests

   !> GAMESS's Bohr, matching the emitter.
   real(dp), parameter :: ANG = 1.0_dp/0.52917724924_dp

   !> The file carries ten decimals, so this is what a round trip can preserve.
   real(dp), parameter :: FORMAT_TOL = 1.0e-9_dp

contains

   subroutine collect_mqc_efp_read_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("efp_round_trip", test_round_trip), &
                  new_unittest("efp_net_charge", test_net_charge), &
                  new_unittest("efp_polarizabilities", test_polarizabilities), &
                  new_unittest("efp_projection_basis", test_projection_basis), &
                  new_unittest("efp_missing_file", test_missing_file) &
                  ]
   end subroutine collect_mqc_efp_read_tests

   subroutine water(pot, err)
      !! One water potential, in the basis the emitter can write
      type(efp_potential_t), intent(out) :: pot
      type(error_t), intent(inout) :: err

      real(dp) :: c(3, 3)
      integer :: z(3)
      character(len=2) :: symbols(3)

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])
      call make_efp_potential(z, symbols, c, "6-31g*", "WATER", pot, err)
   end subroutine water

   subroutine test_round_trip(error)
      type(error_type), allocatable, intent(out) :: error

      type(efp_potential_t) :: pot
      type(efp_fragment_t) :: frag
      type(error_t) :: err
      character(len=*), parameter :: path = "test_efp_round_trip.efp"
      integer :: i, k

      call water(pot, err)
      call check(error,.not. err%has_error(), "building the potential failed")
      if (allocated(error)) return
      call write_efp_potential(pot, path, err)
      call check(error,.not. err%has_error(), "writing the potential failed")
      if (allocated(error)) return

      call read_efp_potential(path, frag, err)
      call check(error,.not. err%has_error(), "reading it back failed")
      if (allocated(error)) return

      call check(error, frag%n_points == pot%n_points, "point count changed")
      if (allocated(error)) return
      call check(error, frag%n_atoms == pot%n_atoms, "atom count changed")
      if (allocated(error)) return
      call check(error, frag%multiplicity == pot%multiplicity, "multiplicity changed")
      if (allocated(error)) return

      do i = 1, frag%n_points
         do k = 1, 3
            call check(error, frag%points(k, i), pot%points(k, i), thr=FORMAT_TOL, &
                       message="an expansion point moved")
            if (allocated(error)) return
         end do
         call check(error, frag%q_elec(i), pot%q_elec(i), thr=FORMAT_TOL, &
                    message="an electronic monopole changed")
         if (allocated(error)) return
         call check(error, frag%q_nuc(i), pot%q_nuc(i), thr=FORMAT_TOL, &
                    message="a nuclear monopole changed")
         if (allocated(error)) return
         do k = 1, 3
            call check(error, frag%dipole(k, i), pot%dipole(k, i), thr=FORMAT_TOL, &
                       message="a dipole component changed")
            if (allocated(error)) return
         end do
         do k = 1, 6
            call check(error, frag%quadrupole(k, i), pot%quadrupole(k, i), &
                       thr=FORMAT_TOL, message="a quadrupole component changed")
            if (allocated(error)) return
         end do
         do k = 1, 10
            call check(error, frag%octopole(k, i), pot%octopole(k, i), &
                       thr=FORMAT_TOL, message="an octupole component changed")
            if (allocated(error)) return
         end do
      end do

      ! Both damping fits are written, so both must come back. A potential read
      ! without them would silently lose the charge-penetration term rather than
      ! fail, which is the kind of thing that shows up as a small energy error.
      call check(error, frag%has_screen, "the Gaussian screening fit was not read")
      if (allocated(error)) return
      call check(error, frag%has_screen2, "the exponential screening fit was not read")
      if (allocated(error)) return
      do i = 1, frag%n_points
         call check(error, frag%screen(i), pot%screen(i), thr=FORMAT_TOL, &
                    message="a Gaussian screening exponent changed")
         if (allocated(error)) return
         call check(error, frag%screen2(i), pot%screen2(i), thr=FORMAT_TOL, &
                    message="an exponential screening exponent changed")
         if (allocated(error)) return
      end do

      call frag%destroy()
      call pot%destroy()
      call delete(path)
   end subroutine test_round_trip

   subroutine test_net_charge(error)
      !! The monopoles must sum to the charge the fragment was computed for
      !!
      !! Neutral water here, so they must sum to zero. This is the cheapest check
      !! that records were split correctly: a monopole record carries two numbers
      !! and reading the wrong one, or reading them in the wrong order, leaves a
      !! fragment with a net charge of eight.
      type(error_type), allocatable, intent(out) :: error

      type(efp_potential_t) :: pot
      type(efp_fragment_t) :: frag
      type(error_t) :: err
      character(len=*), parameter :: path = "test_efp_charge.efp"

      call water(pot, err)
      call check(error,.not. err%has_error(), "building the potential failed")
      if (allocated(error)) return
      call write_efp_potential(pot, path, err)
      call check(error,.not. err%has_error(), "writing the potential failed")
      if (allocated(error)) return
      call read_efp_potential(path, frag, err)
      call check(error,.not. err%has_error(), "reading it back failed")
      if (allocated(error)) return

      call check(error, frag%net_charge(), 0.0_dp, thr=1.0e-8_dp, &
                 message="a neutral fragment did not read back neutral")
      if (allocated(error)) return

      call frag%destroy()
      call pot%destroy()
      call delete(path)
   end subroutine test_net_charge

   subroutine test_polarizabilities(error)
      !! The static and dynamic polarizabilities, checked against each other
      !!
      !! The lowest tabulated frequency is 0.0028 a.u., which is nearly zero, so the
      !! dynamic tensor there must nearly equal the static one. That is the check
      !! that both sections were read into the same slots -- and they are read by
      !! separate code, because their records differ: the static labels are one token
      !! and the dynamic two, and the nine tensor numbers are not row-major but
      !! GAMESS's own slot order, diagonal first with the off-diagonals transposed.
      !!
      !! Reading either wrongly gives a tensor whose trace is negative, so the trace
      !! is asserted positive as well: a polarizability cannot be.
      type(error_type), allocatable, intent(out) :: error

      type(efp_potential_t) :: pot
      type(efp_fragment_t) :: frag
      type(error_t) :: err
      character(len=*), parameter :: path = "test_efp_pol.efp"
      real(dp) :: static_iso, dynamic_iso
      integer :: i

      call water(pot, err)
      call check(error,.not. err%has_error(), "building the potential failed")
      if (allocated(error)) return
      call write_efp_potential(pot, path, err)
      call read_efp_potential(path, frag, err)
      call check(error,.not. err%has_error(), "reading the potential failed")
      if (allocated(error)) return

      call check(error, frag%has_static_pol, "the static polarizabilities were not read")
      if (allocated(error)) return
      call check(error, frag%has_dynamic, "the dynamic polarizabilities were not read")
      if (allocated(error)) return
      call check(error, frag%n_pol == frag%n_lmo, &
                 "the two polarizability sections disagree on the orbital count")
      if (allocated(error)) return
      call check(error, frag%n_freq == 12, "expected twelve frequencies")
      if (allocated(error)) return

      ! The two higher dispersion tensor sets, read flat. Only their extents are
      ! asserted: the slot order is not established, so there is nothing else here
      ! that could be checked without assuming the answer.
      call check(error, frag%has_dipquad, "the dipole-quadrupole tensors were not read")
      if (allocated(error)) return
      call check(error, frag%n_dipquad == 27, "expected 27 dipole-quadrupole values")
      if (allocated(error)) return
      call check(error, frag%has_quadquad, "the quadrupole-quadrupole tensors were not read")
      if (allocated(error)) return
      call check(error, frag%n_quadquad == 81, "expected 81 quadrupole-quadrupole values")
      if (allocated(error)) return
      call check(error, size(frag%dipquad, 2) == frag%n_lmo .and. &
                 size(frag%dipquad, 3) == frag%n_freq, &
                 "the dipole-quadrupole extents disagree with the dipole section")
      if (allocated(error)) return
      ! Not all zero: a section that parsed into nothing would still have the right
      ! shape, and that is the failure this catches.
      call check(error, maxval(abs(frag%quadquad)) > 1.0e-6_dp, &
                 "the quadrupole-quadrupole tensors read as all zero")
      if (allocated(error)) return

      do i = 1, frag%n_pol
         static_iso = (frag%static_pol(1, 1, i) + frag%static_pol(2, 2, i) &
                       + frag%static_pol(3, 3, i))/3.0_dp
         dynamic_iso = (frag%dyn_pol(1, 1, i, 1) + frag%dyn_pol(2, 2, i, 1) &
                        + frag%dyn_pol(3, 3, i, 1))/3.0_dp
         call check(error, static_iso > 0.0_dp, &
                    "a static polarizability came back with a negative trace")
         if (allocated(error)) return
         call check(error, dynamic_iso, static_iso, thr=1.0e-3_dp, &
                    message="the dynamic polarizability at the lowest frequency "// &
                    "does not match the static one")
         if (allocated(error)) return
         ! And the centroids the two sections carry are the same points.
         call check(error, maxval(abs(frag%pol_points(:, i) - frag%centroids(:, i))) &
                    < 1.0e-9_dp, "the two sections disagree on a centroid")
         if (allocated(error)) return
      end do

      call frag%destroy()
      call pot%destroy()
      call delete(path)
   end subroutine test_polarizabilities

   subroutine test_projection_basis(error)
      !! The projection basis, recovered against the basis it was written from
      !!
      !! GAMESS folds primitive normalization into the coefficients it prints, so the
      !! reader divides `gamess_primitive_norm` back out. This checks that the result
      !! is the *original* 6-31G* contraction -- exponents and coefficients both --
      !! by comparing against the same JSON the writer read.
      !!
      !! That is the test that matters for the inter-fragment overlaps to come: a
      !! basis rebuilt with the normalization still folded in would produce overlaps
      !! that are wrong by a smooth factor per primitive, which is exactly the kind
      !! of error that yields plausible energies.
      !!
      !! Shell counts differ by construction: an `L` in the file is one shell there
      !! and two here, matching how the basis file itself holds a shared-exponent sp
      !! pair, so the two sets line up shell for shell.
      type(error_type), allocatable, intent(out) :: error

      type(efp_potential_t) :: pot
      type(efp_fragment_t) :: frag
      type(error_t) :: err
      type(molecular_basis_type) :: basis
      character(len=:), allocatable :: path_json
      character(len=*), parameter :: path = "test_efp_basis.efp"
      character(len=2) :: symbols(3)
      integer :: sh, k, at, want, first
      real(dp) :: worst_e, worst_c

      symbols = ["O ", "H ", "H "]
      call water(pot, err)
      call check(error,.not. err%has_error(), "building the potential failed")
      if (allocated(error)) return
      call write_efp_potential(pot, path, err)
      call read_efp_potential(path, frag, err)
      call check(error,.not. err%has_error(), "reading the potential failed")
      if (allocated(error)) return
      call check(error, frag%has_basis, "the projection basis was not read")
      if (allocated(error)) return

      call find_basis_file("6-31g*", path_json, err)
      call check(error,.not. err%has_error(), "cannot find the basis file")
      if (allocated(error)) return
      call build_molecular_basis_json(path_json, symbols, basis, err)
      call check(error,.not. err%has_error(), "cannot read the basis file")
      if (allocated(error)) return

      want = 0
      do at = 1, 3
         want = want + basis%elements(at)%nshells
      end do
      call check(error, frag%n_shells == want, &
                 "the recovered shell count does not match the basis")
      if (allocated(error)) return

      worst_e = 0.0_dp
      worst_c = 0.0_dp
      sh = 0
      do at = 1, 3
         do k = 1, basis%elements(at)%nshells
            sh = sh + 1
            call check(error, frag%shell_l(sh) == basis%elements(at)%shells(k)%ang_mom, &
                       "a recovered shell has the wrong angular momentum")
            if (allocated(error)) return
            call check(error, frag%shell_nprim(sh) == &
                       size(basis%elements(at)%shells(k)%exponents), &
                       "a recovered shell has the wrong primitive count")
            if (allocated(error)) return
            first = frag%shell_first(sh)
            worst_e = max(worst_e, maxval(abs( &
                                          frag%prim_expo(first:first + frag%shell_nprim(sh) - 1) &
                                          - basis%elements(at)%shells(k)%exponents)))
            worst_c = max(worst_c, maxval(abs( &
                                          frag%prim_coef(first:first + frag%shell_nprim(sh) - 1) &
                                          - basis%elements(at)%shells(k)%coefficients)))
         end do
      end do
      ! The file carries ten decimals on an exponent and eight on a coefficient, so
      ! this is the format's precision rather than the arithmetic's.
      call check(error, worst_e < 1.0e-6_dp, "a recovered exponent does not match")
      if (allocated(error)) return
      call check(error, worst_c < 1.0e-6_dp, "a recovered coefficient does not match")

      call frag%destroy()
      call pot%destroy()
      call delete(path)
   end subroutine test_projection_basis

   subroutine test_missing_file(error)
      !! A path that is not there is an error, not a fragment full of zeros
      type(error_type), allocatable, intent(out) :: error

      type(efp_fragment_t) :: frag
      type(error_t) :: err

      call read_efp_potential("no_such_potential_exists.efp", frag, err)
      call check(error, err%has_error(), "a missing file was read without complaint")
   end subroutine test_missing_file

   subroutine delete(path)
      character(len=*), intent(in) :: path

      integer :: unit, stat

      open (newunit=unit, file=path, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete")
   end subroutine delete

end module test_mqc_efp_read

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_efp_read, only: collect_mqc_efp_read_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_efp_read", collect_mqc_efp_read_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
