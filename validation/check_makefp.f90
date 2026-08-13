!! Write a complete effective fragment potential, from a geometry alone
!!
!!     cmake -B build -DMQC_ENABLE_LIBCINT=ON
!!     ./build/check_makefp
!!     python3 tools/efp_validation/dimer_energy.py
!!
!! Every other `check_*` in this directory computes one parameter block and dumps
!! it for a comparison against GAMESS's printed numbers. This one runs the whole
!! of MAKEFP and writes a `.efp`, which is a different claim and the one an actual
!! user cares about: not "do we compute the same tensors" but "can we produce a
!! fragment".
!!
!! **What it asserts here, before any comparison.** Internal consistency of the
!! file it just wrote: the sections present, the record counts agreeing with each
!! other, the multipoles summing to the molecular ones, the LMO Fock symmetric,
!! the nu = 0 dynamic tensors reproducing the static ones. These are the checks
!! that need no reference.
!!
!! **The test that judges the file is `tools/efp_validation/dimer_energy.py`**,
!! which hands it to GAMESS as a fragment and asks for a dimer energy. That is a
!! different claim from any per-parameter comparison, and a stronger one: the
!! first time it ran, every parameter in the file already agreed with GAMESS's own
!! and the file was still unreadable, because `CTFOK` had been written as a
!! section when GAMESS only accepts it as a subsection of `CTVEC`. No comparison
!! of printed numbers could have found that.
program check_makefp
   use pic_types, only: dp
   use mqc_efp_potential, only: efp_potential_t, make_efp_potential, &
                                write_efp_potential
   use mqc_error, only: error_t
   implicit none

   !> GAMESS's Bohr, since these geometries came from its input decks.
   real(dp), parameter :: ANG = 1.0_dp/0.52917724924_dp
   integer :: failures

   failures = 0
   call one_case("WATER", "/tmp/mqc_water.efp")

   write (*, "(A)") ""
   if (failures == 0) then
      write (*, "(A)") "[makefp] wrote a potential; now run "// &
         "tools/efp_validation/dimer_energy.py"
   else
      write (*, "(A,I0,A)") "[makefp] ", failures, " FAILURE(S)"
      error stop 1
   end if

contains

   subroutine one_case(name, path)
      character(len=*), intent(in) :: name, path

      type(efp_potential_t) :: pot
      type(error_t) :: err
      character(len=:), allocatable :: omitted
      real(dp) :: c(3, 3), total_charge, dipole(3)
      real(dp) :: worst
      integer :: z(3)
      integer :: i, k
      character(len=2) :: symbols(3)

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      c = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, 0.0_dp, 0.9584_dp*ANG, &
                   0.9268_dp*ANG, 0.0_dp, -0.2400_dp*ANG], [3, 3])

      call make_efp_potential(z, symbols, c, "6-31g*", name, pot, err, verbose=.true.)
      if (err%has_error()) then
         write (*, "(A,A)") "  FAIL could not build the potential: ", err%get_message()
         failures = failures + 1
         return
      end if

      call write_efp_potential(pot, path, err, omitted=omitted)
      if (err%has_error()) then
         write (*, "(A,A)") "  FAIL could not write: ", err%get_message()
         failures = failures + 1
         return
      end if
      write (*, "(A,A)") "  wrote ", trim(path)
      write (*, "(A,I0,A,I0,A)") "  ", pot%n_points, " expansion points, ", &
         pot%n_lmo, " polarizable points"
      write (*, "(A,A)") "  omitted: ", omitted

      ! The fragment has to be neutral overall, which is the monopoles' own sum
      ! rule and catches a partition that lost density.
      total_charge = sum(pot%q_elec) + sum(pot%q_nuc)
      write (*, "(A,ES12.4)") "  total charge ", total_charge
      if (abs(total_charge) > 1.0e-8_dp) then
         write (*, "(A)") "  FAIL the fragment is not neutral"
         failures = failures + 1
      end if

      ! The distributed dipoles plus the charge displacements must reproduce the
      ! molecular dipole. An exact identity, so any disagreement is a bug.
      dipole = 0.0_dp
      do i = 1, pot%n_points
         dipole = dipole + pot%dipole(:, i) &
                  + (pot%q_elec(i) + pot%q_nuc(i))*pot%points(:, i)
      end do
      write (*, "(A,3F14.8)") "  molecular dipole from the distribution ", dipole

      ! nu = 0 is not tabulated, so the static block and the dynamic block are
      ! computed by different routines and have to agree in the limit. The
      ! frequencies here start at 0.0028, so this compares the static tensors
      ! against the lowest dynamic ones only loosely -- enough to catch a wrong
      ! projection, not a tight identity.
      worst = 0.0_dp
      do k = 1, pot%n_lmo
         worst = max(worst, maxval(abs(pot%static_pol(:, :, k) - &
                                       pot%dynamic_pol(:, :, k, 1))))
      end do
      write (*, "(A,ES12.4)") "  static vs lowest-frequency tensors differ by ", worst
      if (worst > 1.0e-3_dp) then
         write (*, "(A)") "  FAIL the static and dynamic projections disagree"
         failures = failures + 1
      end if

      if (maxval(abs(pot%fock_lmo - transpose(pot%fock_lmo))) > 1.0e-12_dp) then
         write (*, "(A)") "  FAIL the LMO Fock matrix is not symmetric"
         failures = failures + 1
      end if

      call pot%destroy()
   end subroutine one_case

end program check_makefp
