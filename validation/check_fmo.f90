program check_fmo
   !! Check FMO2 against the supermolecular answer it is approximating
   !!
   !! Three cases, in increasing order of how much they are allowed to be wrong:
   !!
   !!   **water dimer, two fragments.** Not an approximation. With two
   !!   fragments there is nothing outside the dimer, so its embedding field
   !!   vanishes and the monomer terms cancel algebraically, leaving exactly the
   !!   supermolecular SCF. This must agree with full RHF to SCF convergence,
   !!   and if it does not, fragments are being built wrong -- basis
   !!   concatenation, atom ordering, nuclear repulsion. It tests the plumbing
   !!   with none of the physics.
   !!
   !!   **two waters far apart.** Also two fragments and so also exact, but it
   !!   additionally says the monomer SCC loop settles and that a fragment with
   !!   a negligible neighbour is left alone.
   !!
   !!   **water trimer.** Here FMO2 is a real approximation: what it misses is
   !!   the three-body term and whatever the point-charge embedding fails to
   !!   capture. Published FMO2/RHF errors on water clusters are a few
   !!   millihartree, so the tolerance is set there -- loose enough to be a
   !!   statement about the method rather than about this geometry, tight enough
   !!   that a sign error or a double-counted term could not pass.
   !!
   !! The trimer also reports what an unembedded expansion gives on the same
   !! geometry, because the entire argument for FMO over a plain many-body
   !! expansion is that the embedding earns its cost.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger, info_level
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_fmo, only: fmo_options_t, fmo_result_t, run_fmo2
   implicit none

   integer, parameter :: N_DIM = 3
   !> Water plus ammonia, for the unequal-fragment case
   integer, parameter :: N_MIXED = 7
   !> Atom counts for the bond-cutting partition cases
   integer, parameter :: N_ETHANE = 8
   integer, parameter :: N_RING = 9
   integer, parameter :: N_DIMER = 6
   real(dp), parameter :: ANGSTROM_TO_BOHR = 1.8897261254578281_dp
   integer :: n_bad

   n_bad = 0
   call logger%configure(info_level)

   call dimer_case("water dimer", 0.0_dp, 1.0e-7_dp, n_bad)
   call dimer_case("two waters 12 A apart", 12.0_dp, 1.0e-7_dp, n_bad)
   call mixed_case(1.0e-7_dp, n_bad)
   call separated_case(1.0e-9_dp, n_bad)
   call trimer_case(5.0e-3_dp, n_bad)
   call partition_cases(n_bad)
   call level_ladder(n_bad)

   if (n_bad == 0) then
      call logger%info("")
      call logger%info("FMO2 is exact where it must be and close where it may approximate")
   else
      call logger%error("FMO checks failed")
      stop 1
   end if

contains

   subroutine water_geometry(n_waters, push, z, symbols, coords, owner)
      !! `n_waters` copies of the same water, offset along z by `push` Angstrom
      !! beyond the hydrogen-bonded arrangement
      integer, intent(in) :: n_waters
      real(dp), intent(in) :: push
      integer, allocatable, intent(out) :: z(:), owner(:)
      character(len=2), allocatable, intent(out) :: symbols(:)
      real(dp), allocatable, intent(out) :: coords(:, :)

      real(dp) :: monomer(N_DIM, 3)
      real(dp) :: shift
      integer :: w, k, at

      ! One water, Angstrom.
      monomer = reshape([0.0000_dp, 0.0000_dp, 0.0000_dp, &
                         0.0000_dp, -0.7572_dp, 0.5865_dp, &
                         0.0000_dp, 0.7572_dp, 0.5865_dp], [N_DIM, 3])

      allocate (z(3*n_waters), symbols(3*n_waters), owner(3*n_waters))
      allocate (coords(N_DIM, 3*n_waters))

      at = 0
      do w = 1, n_waters
         ! Roughly hydrogen-bond spacing, plus whatever the caller pushes them
         ! apart by. Not an optimised cluster -- the point is a field that is
         ! definitely there, not a physically interesting one.
         shift = real(w - 1, dp)*(2.9_dp + push)
         do k = 1, 3
            at = at + 1
            z(at) = merge(8, 1, k == 1)
            symbols(at) = merge("O ", "H ", k == 1)
            coords(:, at) = monomer(:, k)
            coords(3, at) = coords(3, at) + shift
            owner(at) = w
         end do
      end do

      coords = coords*ANGSTROM_TO_BOHR
   end subroutine water_geometry

   subroutine dimer_case(label, push, tol, n_bad)
      !! Two fragments, where FMO2 is not an approximation at all
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: push, tol
      integer, intent(inout) :: n_bad

      type(fmo_options_t) :: opts
      type(fmo_result_t) :: fmo, eembe
      type(error_t) :: error
      integer, allocatable :: z(:), owner(:)
      character(len=2), allocatable :: symbols(:)
      real(dp), allocatable :: coords(:, :)
      real(dp) :: exact
      character(len=200) :: line

      call logger%info("")
      call logger%info("== "//label)

      call water_geometry(2, push, z, symbols, coords, owner)
      call supermolecule(z, symbols, coords, opts%basis, exact, n_bad, error)
      if (error%has_error()) return

      opts%basis = "6-31g"
      call run_fmo2(z, symbols, coords, owner, opts, fmo, error)
      if (failed(error, "fmo", n_bad)) return

      write (line, "(a,f18.10)") "   supermolecular RHF   ", exact
      call logger%info(trim(line))
      write (line, "(a,f18.10,a,i0,a)") "   FMO2                 ", fmo%energy, &
         "   (", fmo%outer_iterations, " outer passes)"
      call logger%info(trim(line))

      ! With two fragments the dimer has no external field, so the response
      ! term has nothing to respond to. A nonzero value here means the
      ! embedding is leaking into a calculation that should have none.
      call report("no external field on the dimer", abs(fmo%response_sum), 1.0e-12_dp, n_bad)
      call report("FMO2 reproduces the supermolecule", abs(fmo%energy - exact), tol, n_bad)

      ! The same argument holds for the other expansion, and for the same
      ! reason: whatever the monomers were computed in, they cancel. So this
      ! says the many-body assembly is right too, not just FMO's.
      opts%esp = "ptc"
      opts%expansion = "mbe"
      call run_fmo2(z, symbols, coords, owner, opts, eembe, error)
      if (failed(error, "ee-mbe", n_bad)) return
      write (line, "(a,f18.10)") "   EE-MBE               ", eembe%energy
      call logger%info(trim(line))
      call report("EE-MBE reproduces it too", abs(eembe%energy - exact), tol, n_bad)
   end subroutine dimer_case

   subroutine mixed_case(tol, n_bad)
      !! Two fragments of different sizes, which is where an index bug hides
      !!
      !! Water and ammonia have different numbers of basis functions, so the
      !! dimer's two monomer blocks are unequal and the density box has to grow
      !! to fit the larger fragment. Both of those are code paths a cluster of
      !! identical waters never touches. Still two fragments, so still exact --
      !! which means any bookkeeping slip shows up as an outright disagreement
      !! rather than a plausible small error.
      real(dp), intent(in) :: tol
      integer, intent(inout) :: n_bad

      type(fmo_options_t) :: opts
      type(fmo_result_t) :: fmo
      type(error_t) :: error
      integer :: z(N_MIXED), owner(N_MIXED)
      character(len=2) :: symbols(N_MIXED)
      real(dp) :: coords(N_DIM, N_MIXED)
      real(dp) :: exact
      character(len=200) :: line

      call logger%info("")
      call logger%info("== water and ammonia")

      ! Water, then ammonia displaced along z. Angstrom.
      z = [8, 1, 1, 7, 1, 1, 1]
      symbols = [character(len=2) :: "O", "H", "H", "N", "H", "H", "H"]
      owner = [1, 1, 1, 2, 2, 2, 2]
      coords(:, 1) = [0.0000_dp, 0.0000_dp, 0.0000_dp]
      coords(:, 2) = [0.0000_dp, -0.7572_dp, 0.5865_dp]
      coords(:, 3) = [0.0000_dp, 0.7572_dp, 0.5865_dp]
      coords(:, 4) = [0.0000_dp, 0.0000_dp, 3.0000_dp]
      coords(:, 5) = [0.0000_dp, 0.9377_dp, 3.3816_dp]
      coords(:, 6) = [0.8121_dp, -0.4689_dp, 3.3816_dp]
      coords(:, 7) = [-0.8121_dp, -0.4689_dp, 3.3816_dp]
      coords = coords*ANGSTROM_TO_BOHR

      call supermolecule(z, symbols, coords, "6-31g", exact, n_bad, error)
      if (error%has_error()) return

      opts%basis = "6-31g"
      call run_fmo2(z, symbols, coords, owner, opts, fmo, error)
      if (failed(error, "fmo", n_bad)) return

      write (line, "(a,f18.10)") "   supermolecular RHF   ", exact
      call logger%info(trim(line))
      write (line, "(a,f18.10)") "   FMO2                 ", fmo%energy
      call logger%info(trim(line))
      call report("unequal fragments still reproduce the supermolecule", &
                  abs(fmo%energy - exact), tol, n_bad)
   end subroutine mixed_case

   subroutine separated_case(tol, n_bad)
      !! Three fragments far enough apart that only the embedding is left
      !!
      !! The sharpest test of the embedding operator there is. At 9 Angstrom
      !! the fragments still see each other electrostatically -- an unembedded
      !! expansion is out by 1e-06 here -- but exchange and charge transfer
      !! between them have died, so the three-body term FMO2 cannot represent
      !! has died with them. What is left is the embedding, alone.
      !!
      !! An exact ESP therefore has to reproduce the supermolecule to
      !! essentially SCF precision, and it does: 1e-13. Nothing but a correct
      !! embedding operator gets that. Point charges plateau five orders of
      !! magnitude short however far apart the fragments are moved, because the
      !! residual is the charge approximation itself and not a distance effect.
      real(dp), intent(in) :: tol
      integer, intent(inout) :: n_bad

      type(fmo_options_t) :: opts
      type(fmo_result_t) :: exact_esp, ptc, bare, cutoff, ignored
      type(error_t) :: error
      integer, allocatable :: z(:), owner(:)
      character(len=2), allocatable :: symbols(:)
      real(dp), allocatable :: coords(:, :)
      real(dp) :: exact
      character(len=200) :: line

      call logger%info("")
      call logger%info("== three waters 9 A apart")

      ! water_geometry spaces them at 2.9 + push.
      call water_geometry(3, 6.1_dp, z, symbols, coords, owner)
      call supermolecule(z, symbols, coords, "6-31g", exact, n_bad, error)
      if (error%has_error()) return

      opts%basis = "6-31g"
      opts%esp = "exact"
      ! 9 Angstrom over the summed van der Waals radii is a unitless separation
      ! near 3, so under the default RESPPC of 2 these fragments are distant and
      ! would be approximated. This case is about the exact operator, so it asks
      ! for it: a negative cutoff turns the approximation off.
      opts%resppc = -1.0_dp
      call run_fmo2(z, symbols, coords, owner, opts, exact_esp, error)
      if (failed(error, "fmo exact esp", n_bad)) return

      opts%esp = "ptc"
      opts%far_field = "mulliken"
      call run_fmo2(z, symbols, coords, owner, opts, ptc, error)
      if (failed(error, "fmo ptc", n_bad)) return

      ! And the default cutoff must agree with the point-charge run here, since
      ! at this separation it should be choosing point charges for everything.
      opts%esp = "exact"
      opts%resppc = 2.0_dp
      call run_fmo2(z, symbols, coords, owner, opts, cutoff, error)
      if (failed(error, "fmo default cutoff", n_bad)) return

      opts%esp = "none"
      opts%expansion = "mbe"
      call run_fmo2(z, symbols, coords, owner, opts, bare, error)
      if (failed(error, "plain mbe", n_bad)) return

      ! Every fragment distant, and distant fragments ignored: there is then no
      ! embedding left anywhere, so this has to be plain MBE exactly. It reaches
      ! that through the near/far machinery rather than the no-embedding branch,
      ! which is why it is worth asking.
      opts%esp = "ptc"
      opts%expansion = "mbe"
      opts%far_field = "ignore"
      call run_fmo2(z, symbols, coords, owner, opts, ignored, error)
      if (failed(error, "far field ignored", n_bad)) return
      opts%far_field = "mulliken"

      write (line, "(a,es11.3)") "   exact ESP      err ", exact_esp%energy - exact
      call logger%info(trim(line))
      write (line, "(a,es11.3)") "   ptc, mulliken  err ", ptc%energy - exact
      call logger%info(trim(line))
      write (line, "(a,es11.3)") "   no embedding   err ", bare%energy - exact
      call logger%info(trim(line))

      write (line, "(a,es11.3)") "   default cutoff err ", cutoff%energy - exact
      call logger%info(trim(line))

      call report("exact ESP is exact once 3-body effects have died", &
                  abs(exact_esp%energy - exact), tol, n_bad)

      ! At this separation the cutoff should have classified every fragment as
      ! distant, so it must land exactly on the point-charge answer -- reached
      ! by a different route, through the near/far split rather than through the
      ! all-point-charge branch.
      call report("the default cutoff picks point charges out here", &
                  abs(cutoff%energy - ptc%energy), 1.0e-12_dp, n_bad)
      call report("ignoring every distant fragment leaves plain MBE", &
                  abs(ignored%energy - bare%energy), 1.0e-12_dp, n_bad)

      ! The point of separating them: the field is still real at this distance,
      ! so an expansion that ignores it is visibly wrong. Were that to stop
      ! holding, the case would have gone slack -- everything would agree
      ! because there was nothing left to get right -- and it would no longer
      ! be testing the embedding at all.
      if (abs(bare%energy - exact) > 1.0e-7_dp) then
         write (line, "(a,es9.2,a)") "   ok   the field is still real: ignoring it costs ", &
            abs(bare%energy - exact), " Hartree"
         call logger%info(trim(line))
      else
         call logger%error("   FAIL the fragments are now too far apart to be testing "// &
                           "anything; move them closer")
         n_bad = n_bad + 1
      end if

      ! Point charges cannot get here however far apart the fragments go.
      if (abs(ptc%energy - exact) <= abs(exact_esp%energy - exact)) then
         call logger%error("   FAIL point charges matched the exact ESP at long range, "// &
                           "where the exact one should be exact")
         n_bad = n_bad + 1
      else
         write (line, "(a,es9.2,a)") "   ok   point charges plateau at ", &
            abs(ptc%energy - exact), ", the charge approximation itself"
         call logger%info(trim(line))
      end if
   end subroutine separated_case

   subroutine trimer_case(tol, n_bad)
      !! Three fragments, where the embedding starts to matter
      real(dp), intent(in) :: tol
      integer, intent(inout) :: n_bad

      type(fmo_options_t) :: opts
      type(fmo_result_t) :: embedded, bare, fitted, eembe
      type(error_t) :: error
      integer, allocatable :: z(:), owner(:)
      character(len=2), allocatable :: symbols(:)
      real(dp), allocatable :: coords(:, :)
      real(dp) :: exact
      character(len=200) :: line

      call logger%info("")
      call logger%info("== water trimer")

      call water_geometry(3, 0.0_dp, z, symbols, coords, owner)
      call supermolecule(z, symbols, coords, "6-31g", exact, n_bad, error)
      if (error%has_error()) return

      opts%basis = "6-31g"
      opts%esp = "exact"
      call run_fmo2(z, symbols, coords, owner, opts, embedded, error)
      if (failed(error, "fmo exact esp", n_bad)) return

      opts%esp = "ptc"
      opts%far_field = "chelpg"
      call run_fmo2(z, symbols, coords, owner, opts, fitted, error)
      if (failed(error, "fmo ptc chelpg", n_bad)) return

      opts%esp = "ptc"
      opts%expansion = "mbe"
      opts%far_field = "mulliken"
      call run_fmo2(z, symbols, coords, owner, opts, eembe, error)
      if (failed(error, "ee-mbe", n_bad)) return

      opts%esp = "none"
      opts%expansion = "mbe"
      call run_fmo2(z, symbols, coords, owner, opts, bare, error)
      if (failed(error, "plain mbe", n_bad)) return

      write (line, "(a,f18.10)") "   supermolecular RHF        ", exact
      call logger%info(trim(line))
      write (line, "(a,f18.10,a,es10.2,a,i0,a)") "   FMO2, exact ESP           ", &
         embedded%energy, "   err ", embedded%energy - exact, &
         "   (", embedded%outer_iterations, " outer)"
      call logger%info(trim(line))
      write (line, "(a,f18.10,a,es10.2,a,i0,a)") "   FMO2, ptc chelpg          ", &
         fitted%energy, "   err ", fitted%energy - exact, &
         "   (", fitted%outer_iterations, " outer)"
      call logger%info(trim(line))
      write (line, "(a,f18.10,a,es10.2,a,i0,a)") "   EE-MBE, mulliken          ", &
         eembe%energy, "   err ", eembe%energy - exact, &
         "   (", eembe%outer_iterations, " outer)"
      call logger%info(trim(line))
      write (line, "(a,f18.10,a,es10.2)") "   plain MBE, no embedding   ", &
         bare%energy, "   err ", bare%energy - exact
      call logger%info(trim(line))
      write (line, "(a,f14.8)") "   the response term contributes ", embedded%response_sum
      call logger%info(trim(line))

      call report("exact-ESP FMO2 is close to exact", abs(embedded%energy - exact), tol, n_bad)
      call report("ptc-chelpg FMO2 is close to exact", abs(fitted%energy - exact), &
                  tol, n_bad)

      ! The point of the whole exercise. If embedding did not help, there would
      ! be no reason to prefer this to a plain many-body expansion, and the
      ! monomer SCC loop would be cost for nothing.
      if (abs(embedded%energy - exact) < abs(bare%energy - exact)) then
         write (line, "(a,f6.2,a)") "   ok   embedding helps, by ", &
            abs(bare%energy - exact)/abs(embedded%energy - exact), "x"
         call logger%info(trim(line))
      else
         call logger%error("   FAIL embedding did not beat no embedding, which is the "// &
                           "whole argument for FMO over an unembedded expansion")
         n_bad = n_bad + 1
      end if
   end subroutine trimer_case

   subroutine level_ladder(n_bad)
      !! Climbing the expansion level has to converge on the exact answer
      !!
      !! FMO_n truncates the many-body expansion at n fragments at a time. When
      !! n reaches the number of fragments there is nothing left to truncate:
      !! the corrections telescope and the result is the supermolecular energy,
      !! by the same argument that makes FMO2 exact on two fragments. So the
      !! ladder ends on an equality, not on a tolerance, and every rung below it
      !! should be closer than the one before.
      !!
      !! That upper rung is the real test of the recursion. A wrong many-body
      !! coefficient still gives plausible numbers at level 2, where the formula
      !! is short enough to be right by accident; it cannot survive being asked
      !! to cancel exactly.
      integer, intent(inout) :: n_bad

      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      type(error_t) :: error
      integer, allocatable :: z(:), owner(:)
      character(len=2), allocatable :: symbols(:)
      real(dp), allocatable :: coords(:, :)
      real(dp) :: exact, err, previous
      character(len=200) :: line
      integer :: n_waters, lvl

      do n_waters = 3, 4
         call logger%info("")
         write (line, "(a,i0,a)") "== ", n_waters, " waters, level by level"
         call logger%info(trim(line))

         call water_geometry(n_waters, 0.0_dp, z, symbols, coords, owner)
         call supermolecule(z, symbols, coords, "sto-3g", exact, n_bad, error)
         if (error%has_error()) return

         previous = huge(1.0_dp)
         do lvl = 2, n_waters
            opts%basis = "sto-3g"
            opts%esp = "exact"
            opts%expansion = "fmo"
            opts%resppc = -1.0_dp     !! no cutoff: this is about the expansion
            opts%level = lvl
            call run_fmo2(z, symbols, coords, owner, opts, res, error)
            if (failed(error, "fmo level", n_bad)) return

            err = abs(res%energy - exact)
            write (line, "(a,i0,a,f18.10,a,es11.3)") "   level ", lvl, "   E = ", &
               res%energy, "   err ", err
            call logger%info(trim(line))

            if (err > previous) then
               call logger%error("   FAIL raising the level made it worse")
               n_bad = n_bad + 1
            end if
            previous = err
         end do

         ! The last rung used every fragment at once, so it is not an
         ! approximation and the tolerance is SCF convergence, not chemistry.
         call report("the full expansion is the supermolecule", err, 1.0e-9_dp, n_bad)
      end do
   end subroutine level_ladder

   subroutine partition_cases(n_bad)
      !! A partition that cuts a covalent bond has to be refused, not answered
      !!
      !! There is no capping here -- no adjusted fragment orbitals, no hybrid
      !! orbital projection, no hydrogen caps -- so a fragment with a dangling
      !! valence is a question the method cannot be asked. The reason to test it
      !! rather than trust it is that the failure used to be quiet in exactly
      !! the cases that matter:
      !!
      !!   * Cutting one bond leaves both fragments with an odd electron count,
      !!     and the closed-shell check caught that by accident.
      !!   * Cutting an even number per fragment -- a ring, a double bond --
      !!     left every count even, nothing objected, and cyclopropane split
      !!     into three CH2 came back 0.28 Hartree low. That is 176 kcal/mol
      !!     wearing the shape of an answer.
      !!
      !! The third case is the one that keeps the check honest in the other
      !! direction: a hydrogen bond is not a covalent one, and a tightly bound
      !! water dimer must still be accepted. A check that refused those would be
      !! useless in precisely the systems this method is for.
      integer, intent(inout) :: n_bad

      real(dp) :: ethane(N_DIM, N_ETHANE)
      real(dp) :: ring(N_DIM, N_RING)
      real(dp) :: dimer(N_DIM, N_DIMER)

      call logger%info("")
      call logger%info("== partitions that cut bonds")

      ethane = reshape([0.000_dp, 0.000_dp, 0.768_dp, &
                        -1.019_dp, 0.000_dp, 1.157_dp, &
                        0.510_dp, 0.883_dp, 1.157_dp, &
                        0.510_dp, -0.883_dp, 1.157_dp, &
                        0.000_dp, 0.000_dp, -0.768_dp, &
                        1.019_dp, 0.000_dp, -1.157_dp, &
                        -0.510_dp, -0.883_dp, -1.157_dp, &
                        -0.510_dp, 0.883_dp, -1.157_dp], [N_DIM, N_ETHANE])
      call must_refuse("ethane cut at C-C", [6, 1, 1, 1, 6, 1, 1, 1], ethane, &
                       [1, 1, 1, 1, 2, 2, 2, 2], n_bad)

      ! Three CH2, eight electrons each: every count even, so nothing but a
      ! connectivity check stands between this and a confident wrong answer.
      ring = reshape([0.8718_dp, 0.0000_dp, 0.0000_dp, &
                      1.3818_dp, 0.0000_dp, 0.9500_dp, &
                      1.3818_dp, 0.0000_dp, -0.9500_dp, &
                      -0.4359_dp, 0.7550_dp, 0.0000_dp, &
                      -0.6909_dp, 1.1967_dp, 0.9500_dp, &
                      -0.6909_dp, 1.1967_dp, -0.9500_dp, &
                      -0.4359_dp, -0.7550_dp, 0.0000_dp, &
                      -0.6909_dp, -1.1967_dp, 0.9500_dp, &
                      -0.6909_dp, -1.1967_dp, -0.9500_dp], [N_DIM, N_RING])
      call must_refuse("cyclopropane cut into three CH2", [6, 1, 1, 6, 1, 1, 6, 1, 1], &
                       ring, [1, 1, 1, 2, 2, 2, 3, 3, 3], n_bad)

      ! And the other direction: 2.5 Angstrom is tighter than any hydrogen bond
      ! these tests use, and must still be allowed through.
      dimer = reshape([0.0000_dp, 0.0000_dp, 0.0000_dp, &
                       0.0000_dp, -0.7572_dp, 0.5865_dp, &
                       0.0000_dp, 0.7572_dp, 0.5865_dp, &
                       0.0000_dp, 0.0000_dp, 2.5000_dp, &
                       0.0000_dp, -0.7572_dp, 3.0865_dp, &
                       0.0000_dp, 0.7572_dp, 3.0865_dp], [N_DIM, N_DIMER])
      call must_allow("two waters 2.5 A apart", [8, 1, 1, 8, 1, 1], dimer, &
                      [1, 1, 1, 2, 2, 2], n_bad)
   end subroutine partition_cases

   subroutine must_refuse(label, z, coords_ang, owner, n_bad)
      character(len=*), intent(in) :: label
      integer, intent(in) :: z(:), owner(:)
      real(dp), intent(in) :: coords_ang(:, :)
      integer, intent(inout) :: n_bad

      type(error_t) :: error

      call attempt(z, coords_ang, owner, error)
      if (error%has_error()) then
         call logger%info("   ok   "//label//" is refused")
         call error%clear()
      else
         call logger%error("   FAIL "//label//" was answered rather than refused")
         n_bad = n_bad + 1
      end if
   end subroutine must_refuse

   subroutine must_allow(label, z, coords_ang, owner, n_bad)
      character(len=*), intent(in) :: label
      integer, intent(in) :: z(:), owner(:)
      real(dp), intent(in) :: coords_ang(:, :)
      integer, intent(inout) :: n_bad

      type(error_t) :: error

      call attempt(z, coords_ang, owner, error)
      if (error%has_error()) then
         call logger%error("   FAIL "//label//" was refused: "//error%get_message())
         call error%clear()
         n_bad = n_bad + 1
      else
         call logger%info("   ok   "//label//" is allowed, a hydrogen bond is not a cut")
      end if
   end subroutine must_allow

   subroutine attempt(z, coords_ang, owner, error)
      integer, intent(in) :: z(:), owner(:)
      real(dp), intent(in) :: coords_ang(:, :)
      type(error_t), intent(inout) :: error

      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      character(len=2), allocatable :: symbols(:)
      real(dp), allocatable :: coords(:, :)
      integer :: i

      allocate (symbols(size(z)))
      do i = 1, size(z)
         select case (z(i))
         case (6)
            symbols(i) = "C "
         case (7)
            symbols(i) = "N "
         case (8)
            symbols(i) = "O "
         case default
            symbols(i) = "H "
         end select
      end do
      coords = coords_ang*ANGSTROM_TO_BOHR

      opts%basis = "sto-3g"
      call run_fmo2(z, symbols, coords, owner, opts, res, error)
   end subroutine attempt

   subroutine supermolecule(z, symbols, coords, basis, energy, n_bad, error)
      !! The answer FMO2 is trying to reproduce
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      real(dp), intent(out) :: energy
      integer, intent(inout) :: n_bad
      type(error_t), intent(inout) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf

      energy = 0.0_dp
      call build_libcint_molecule(z, symbols, coords, basis, mol, error)
      if (failed(error, "supermolecule", n_bad)) return
      call run_libcint_rhf(mol, sum(z), 200, 1.0e-10_dp, 1.0e-8_dp, .false., scf, error)
      if (failed(error, "supermolecular scf", n_bad)) return
      energy = scf%energy
   end subroutine supermolecule

   function failed(error, what, n_bad) result(bad)
      type(error_t), intent(inout) :: error
      character(len=*), intent(in) :: what
      integer, intent(inout) :: n_bad
      logical :: bad

      bad = error%has_error()
      if (bad) then
         call logger%error("   FAIL "//what//": "//error%get_message())
         n_bad = n_bad + 1
      end if
   end function failed

   subroutine report(what, deviation, tol, n_bad)
      character(len=*), intent(in) :: what
      real(dp), intent(in) :: deviation, tol
      integer, intent(inout) :: n_bad

      character(len=200) :: line

      if (deviation <= tol) then
         write (line, "(a,a,a,es10.3)") "   ok   ", what, "   ", deviation
         call logger%info(trim(line))
      else
         write (line, "(a,a,a,es10.3,a,es10.3)") "   FAIL ", what, "   ", deviation, &
            " > ", tol
         call logger%error(trim(line))
         n_bad = n_bad + 1
      end if
   end subroutine report

end program check_fmo
