!! Dropping grid points where the density is negligible
module test_mqc_xc_prune
   !! `xc_prune_small_rho` removes the points where `abs(rho*w)` on a density is
   !! at or below the cutoff divided by the grid size, and nothing else about the quadrature may change.
   !! Checked on a water dimer, PBE, def2-SVP, grid level 3: that points go, that
   !! every per-point array of the grid is compacted together and keeps its order,
   !! that the electrons the grid integrates to move by exactly what was reported
   !! dropped, that a cutoff of zero keeps everything, that a second call on the
   !! same context does nothing, and that an SCF started on a pruned context lands
   !! on the unpruned energy to well under a millihartree.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_xc, only: xc_context_t, xc_context_create, xc_add_potential, &
                         xc_prune_small_rho, xc_available
   implicit none
   private
   public :: collect_mqc_xc_prune

   real(dp), parameter :: ANG = 1.8897261254578281_dp
   integer, parameter :: N_ELEC = 20
   integer, parameter :: N_ATOMS = 6

contains

   subroutine collect_mqc_xc_prune(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("pruning_drops_points_and_keeps_the_grid_consistent", test_prune), &
                  new_unittest("a_zero_cutoff_keeps_every_point", test_zero_cutoff), &
                  new_unittest("the_scf_prunes_once_and_keeps_the_energy", test_scf) &
                  ]
   end subroutine collect_mqc_xc_prune

   subroutine dimer(mol, err)
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: err
      real(dp) :: c(3, N_ATOMS)

      c = reshape([-1.551007_dp, -0.114520_dp, 0.0_dp, &
                   -1.934259_dp, 0.762503_dp, 0.0_dp, &
                   -0.599677_dp, 0.040712_dp, 0.0_dp, &
                   1.350625_dp, 0.111469_dp, 0.0_dp, &
                   1.680398_dp, -0.373741_dp, -0.758561_dp, &
                   1.680398_dp, -0.373741_dp, 0.758561_dp], [3, N_ATOMS])*ANG
      call build_czt_molecule([8, 1, 1, 8, 1, 1], ["O ", "H ", "H ", "O ", "H ", "H "], c, &
                              "def2-svp", mol, err)
   end subroutine dimer

   subroutine converged_density(mol, density, energy, error)
      !! An unpruned PBE SCF, for the density the prune is asked about
      type(czt_molecule_t), intent(in) :: mol
      real(dp), allocatable, intent(out) :: density(:, :)
      real(dp), intent(out) :: energy
      type(error_type), allocatable, intent(out) :: error
      type(xc_context_t) :: xc
      type(rhf_result_t) :: res
      type(error_t) :: err

      call xc_context_create(mol, "pbe", xc, err, level=3)
      call check(error,.not. err%has_error(), "the context must build: "//err%get_full_trace())
      if (allocated(error)) return
      call run_czt_rhf(mol, N_ELEC, 100, 1.0e-10_dp, 1.0e-8_dp, .false., res, err, xc=xc)
      call check(error,.not. err%has_error(), "the SCF must run: "//err%get_full_trace())
      if (allocated(error)) return
      density = res%density
      energy = res%energy
      call xc%destroy()
   end subroutine converged_density

   subroutine test_prune(error)
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(xc_context_t) :: full, cut
      type(error_t) :: err
      real(dp), allocatable :: density(:, :), v(:, :)
      real(dp) :: energy, e_full, n_full, n_cut, e_cut, dropped
      integer :: n0, n1, i, j, maxbad

      if (.not. xc_available()) return
      call dimer(mol, err)
      call check(error,.not. err%has_error(), "the dimer must build: "//err%get_full_trace())
      if (allocated(error)) return
      call converged_density(mol, density, energy, error)
      if (allocated(error)) return

      call xc_context_create(mol, "pbe", full, err, level=3)
      call xc_context_create(mol, "pbe", cut, err, level=3, small_rho_cutoff=1.0e-7_dp)
      call check(error,.not. err%has_error(), "the contexts must build: "//err%get_full_trace())
      if (allocated(error)) return
      n0 = full%grid%n_points
      allocate (v(size(density, 1), size(density, 2)))
      call xc_add_potential(full, mol, density, v, e_full, n_full, err)

      call xc_prune_small_rho(cut, mol, density, N_ELEC, err, dropped)
      call check(error,.not. err%has_error(), "the prune must run: "//err%get_full_trace())
      if (allocated(error)) return
      n1 = cut%grid%n_points
      call check(error, n1 < n0, "pruning must remove points from a water dimer grid")
      if (allocated(error)) return
      call check(error, n1 > n0/2, "pruning must not remove most of the grid")
      if (allocated(error)) return

      ! Every per-point array has the new length.
      call check(error, size(cut%grid%coords, 2) == n1 .and. size(cut%grid%weights) == n1 &
                 .and. size(cut%grid%quad_weights) == n1 .and. size(cut%grid%atom) == n1, &
                 "every per-point array must be compacted to the same length")
      if (allocated(error)) return

      ! The kept points are the originals, in order: walk the old grid once.
      j = 0
      maxbad = 0
      do i = 1, n1
         do
            j = j + 1
            if (j > n0) exit
            if (all(full%grid%coords(:, j) == cut%grid%coords(:, i))) exit
         end do
         if (j > n0) then
            maxbad = maxbad + 1
            exit
         end if
         if (full%grid%weights(j) /= cut%grid%weights(i) .or. &
             full%grid%quad_weights(j) /= cut%grid%quad_weights(i) .or. &
             full%grid%atom(j) /= cut%grid%atom(i)) maxbad = maxbad + 1
      end do
      call check(error, maxbad == 0, "kept points must be the original ones, in order, with "// &
                 "their own weights and owners")
      if (allocated(error)) return

      ! What the grid integrates to moves by what was reported dropped.
      call xc_add_potential(cut, mol, density, v, e_cut, n_cut, err)
      call check(error, abs((n_full - n_cut) - dropped) < 1.0e-9_dp, &
                 "the electron count must fall by exactly the reported dropped amount")
      if (allocated(error)) return
      call check(error, dropped >= 0.0_dp .and. dropped < real(n0 - n1, dp)*1.0e-7_dp/real(n0, dp) + 1.0e-12_dp, &
                 "each dropped point carries at most the cutoff")
      if (allocated(error)) return
      ! An electron in the tail carries well under a hartree of XC energy, so the
      ! energy moves by less than the dropped count in hartree.
      call check(error, abs(e_full - e_cut) < dropped, &
                 "the XC energy must move by less than the dropped electrons in hartree")

      call full%destroy()
      call cut%destroy()
      call mol%destroy()
   end subroutine test_prune

   subroutine test_zero_cutoff(error)
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(xc_context_t) :: xc
      type(error_t) :: err
      real(dp), allocatable :: density(:, :)
      real(dp) :: energy, dropped

      if (.not. xc_available()) return
      call dimer(mol, err)
      call check(error,.not. err%has_error(), "the dimer must build: "//err%get_full_trace())
      if (allocated(error)) return
      call converged_density(mol, density, energy, error)
      if (allocated(error)) return
      call xc_context_create(mol, "pbe", xc, err, level=3, small_rho_cutoff=0.0_dp)
      call xc_prune_small_rho(xc, mol, density, N_ELEC, err, dropped)
      call check(error, dropped == 0.0_dp, "nothing is dropped at a cutoff of zero")
      call xc%destroy()
      call mol%destroy()
   end subroutine test_zero_cutoff

   subroutine test_scf(error)
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(xc_context_t) :: xc0, xc1
      type(rhf_result_t) :: r0, r1
      type(error_t) :: err
      integer :: n_full

      if (.not. xc_available()) return
      call dimer(mol, err)
      call check(error,.not. err%has_error(), "the dimer must build: "//err%get_full_trace())
      if (allocated(error)) return
      call xc_context_create(mol, "pbe", xc0, err, level=3)
      call xc_context_create(mol, "pbe", xc1, err, level=3, small_rho_cutoff=1.0e-7_dp)
      n_full = xc0%grid%n_points
      call run_czt_rhf(mol, N_ELEC, 100, 1.0e-10_dp, 1.0e-8_dp, .false., r0, err, xc=xc0)
      call run_czt_rhf(mol, N_ELEC, 100, 1.0e-10_dp, 1.0e-8_dp, .false., r1, err, xc=xc1)
      call check(error,.not. err%has_error(), "the SCFs must run: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, xc0%grid%n_points == n_full, "a zero-cutoff context keeps its grid")
      if (allocated(error)) return
      call check(error, xc1%grid%n_points < n_full, "the SCF must prune its grid")
      if (allocated(error)) return
      call check(error, xc1%grid_pruned, "the context must record that it was pruned")
      if (allocated(error)) return
      n_full = xc1%grid%n_points
      call xc_prune_small_rho(xc1, mol, r1%density, N_ELEC, err)
      call check(error, xc1%grid%n_points == n_full, "a second prune of the same context does nothing")
      if (allocated(error)) return
      call check(error, abs(r0%energy - r1%energy) < 1.0e-3_dp, &
                 "the pruned energy must agree with the unpruned one to a millihartree")
      call xc0%destroy()
      call xc1%destroy()
      call mol%destroy()
   end subroutine test_scf

end module test_mqc_xc_prune

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_xc_prune, only: collect_mqc_xc_prune
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_xc_prune", collect_mqc_xc_prune)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
