program check_projection
   !! Check the basis-set projection guess against what it has to be
   !!
   !! Three properties, and they are checkable without a reference implementation
   !! because each is exact by construction rather than approximately right:
   !!
   !!   1. **Trace.** `Tr[D S_BB]` is the electron count. This is what the naive
   !!      density projection gets wrong and is the reason the orbitals are
   !!      projected instead.
   !!   2. **Idempotency.** `D S D = 2 D` for a closed shell with `D = 2 C C^T`
   !!      and `C^T S C = I`. A density that fails this is not a density.
   !!   3. **Identity.** Projecting a basis onto *itself* must return the density
   !!      it started from. This is the one that catches a wrong cross-overlap
   !!      block, a mis-shifted `env` pointer or a transposed projection, none of
   !!      which the first two would notice.
   !!
   !! Then the point of the exercise: does starting from the projected density
   !! take fewer iterations than the default guess.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_logger, only: logger => global_logger, info_level
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_projection, only: project_occupied, cross_overlap
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, SCF_GUESS_CORE, SCF_GUESS_SAD
   implicit none

   integer, parameter :: N_DIM = 3
   integer :: n_bad

   n_bad = 0
   call logger%configure(info_level)

   ! Water, which converges easily -- this case is about the projection being
   ! right, not about it being needed.
   call run_case("water", [character(len=2) :: "O", "H", "H"], &
                 reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, -1.4308_dp, 1.1078_dp, &
                          0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                 [8, 1, 1], 10, "sto-3g", "6-31g", n_bad)

   ! Two heavy atoms and a longer basis, so the cross-overlap has off-diagonal
   ! structure a single atom would not exercise.
   call run_case("ethene", [character(len=2) :: "C", "C", "H", "H", "H", "H"], &
                 reshape([0.0_dp, 0.0_dp, 1.2637_dp, &
                          0.0_dp, 0.0_dp, -1.2637_dp, &
                          0.0_dp, 1.7554_dp, 2.3286_dp, &
                          0.0_dp, -1.7554_dp, 2.3286_dp, &
                          0.0_dp, 1.7554_dp, -2.3286_dp, &
                          0.0_dp, -1.7554_dp, -2.3286_dp], [N_DIM, 6]), &
                 [6, 6, 1, 1, 1, 1], 16, "sto-3g", "6-31g", n_bad)

   if (n_bad == 0) then
      write (*, "(a)") ""
      write (*, "(a)") "projection guess is idempotent, normalised, and exact onto itself"
   else
      write (*, "(a,i0,a)") "", n_bad, " check(s) failed"
      stop 1
   end if

contains

   subroutine run_case(label, symbols, coords, charges, nelec, small, target_basis, n_bad)
      character(len=*), intent(in) :: label, small, target_basis
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: charges(:)
      integer, intent(in) :: nelec
      integer, intent(inout) :: n_bad

      type(czt_molecule_t) :: mol_s, mol_t
      type(rhf_result_t) :: scf_s, scf_t_proj, scf_t_core
      type(error_t) :: error
      real(dp), allocatable :: d_proj(:, :), s_t(:, :), s_x(:, :), d_self(:, :)
      real(dp), allocatable :: ds(:, :), dsd(:, :)
      real(dp) :: trace, worst_idem, worst_self
      integer :: n_occ, i

      n_occ = nelec/2
      write (*, "(a)") ""
      write (*, "(a,a,a,a,a,a)") "== ", label, ": ", small, " -> ", target_basis

      call build_czt_molecule(charges, symbols, coords, small, mol_s, error)
      if (failed(error, "small molecule", n_bad)) return
      call build_czt_molecule(charges, symbols, coords, target_basis, mol_t, error)
      if (failed(error, "target molecule", n_bad)) return

      call run_czt_rhf(mol_s, nelec, 100, 1.0e-9_dp, 1.0e-7_dp, .false., scf_s, error, &
                       guess=SCF_GUESS_CORE)
      if (failed(error, "small-basis SCF", n_bad)) return
      write (*, "(a,i0,a,i0,a,f18.10)") "  ", mol_s%nao, " -> ", mol_t%nao, &
         " functions,  small-basis E = ", scf_s%energy

      ! ---- 3. identity: the same basis both sides must reproduce its density ----
      call project_occupied(mol_s, mol_s, scf_s%orbitals, n_occ, d_self, error)
      if (failed(error, "self-projection", n_bad)) return
      worst_self = maxval(abs(d_self - scf_s%density))
      call report("self-projection reproduces its own density", worst_self, 1.0e-10_dp, n_bad)

      ! ---- the projection under test ----
      call project_occupied(mol_t, mol_s, scf_s%orbitals, n_occ, d_proj, error)
      if (failed(error, "projection", n_bad)) return
      call cross_overlap(mol_t, mol_s, s_t, s_x, error)
      if (failed(error, "cross overlap", n_bad)) return

      ! ---- 1. trace ----
      allocate (ds(mol_t%nao, mol_t%nao))
      call pic_gemm(d_proj, s_t, ds)
      trace = 0.0_dp
      do i = 1, mol_t%nao
         trace = trace + ds(i, i)
      end do
      call report("Tr[D S] is the electron count", abs(trace - real(nelec, dp)), &
                  1.0e-10_dp, n_bad)

      ! ---- 2. idempotency: D S D = 2 D ----
      allocate (dsd(mol_t%nao, mol_t%nao))
      call pic_gemm(ds, d_proj, dsd)
      worst_idem = maxval(abs(dsd - 2.0_dp*d_proj))
      call report("D S D = 2 D", worst_idem, 1.0e-10_dp, n_bad)

      ! ---- what it is for ----
      call run_czt_rhf(mol_t, nelec, 100, 1.0e-9_dp, 1.0e-7_dp, .false., scf_t_core, &
                       error, guess=SCF_GUESS_CORE)
      if (failed(error, "target SCF from core", n_bad)) return
      call run_czt_rhf(mol_t, nelec, 100, 1.0e-9_dp, 1.0e-7_dp, .false., scf_t_proj, &
                       error, guess=SCF_GUESS_SAD, guess_density=d_proj)
      if (failed(error, "target SCF from projection", n_bad)) return

      write (*, "(a,i0,a,f18.10)") "  from core guess:        ", scf_t_core%iterations, &
         " iterations, E = ", scf_t_core%energy
      write (*, "(a,i0,a,f18.10)") "  from projected density: ", scf_t_proj%iterations, &
         " iterations, E = ", scf_t_proj%energy
      call report("both starting points reach the same energy", &
                  abs(scf_t_core%energy - scf_t_proj%energy), 1.0e-8_dp, n_bad)
   end subroutine run_case

   function failed(error, what, n_bad) result(bad)
      type(error_t), intent(inout) :: error
      character(len=*), intent(in) :: what
      integer, intent(inout) :: n_bad
      logical :: bad

      bad = error%has_error()
      if (bad) then
         write (*, "(a,a,a,a)") "  FAIL ", what, ": ", error%get_message()
         n_bad = n_bad + 1
      end if
   end function failed

   subroutine report(what, deviation, tol, n_bad)
      character(len=*), intent(in) :: what
      real(dp), intent(in) :: deviation, tol
      integer, intent(inout) :: n_bad

      if (deviation <= tol) then
         write (*, "(a,a,a,es10.3)") "  ok   ", what, "   ", deviation
      else
         write (*, "(a,a,a,es10.3,a,es10.3)") "  FAIL ", what, "   ", deviation, " > ", tol
         n_bad = n_bad + 1
      end if
   end subroutine report

end program check_projection
