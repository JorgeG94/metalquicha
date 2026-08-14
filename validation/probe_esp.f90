program probe_esp
   !! Is the fast embedding really the slow one?
   !!
   !! The exact ESP needs, for fragment X, the Coulomb operator of every *other*
   !! fragment's density. Built the obvious way that is one supersystem Coulomb
   !! matrix per fragment per outer pass, from a density with X's block zeroed.
   !! [[mqc_libcint_fmo]] instead builds one supersystem matrix from the total
   !! density and subtracts a small one over X's own basis, which is the same
   !! thing only if a four-index integral over X's functions is independent of
   !! which molecule object computed it.
   !!
   !! That identity is the load-bearing assumption of the whole embedding, so it
   !! is checked here directly rather than inferred from an energy coming out
   !! plausible.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger, info_level
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_direct, only: schwarz_bounds, build_fock_direct, direct_stats_t
   use mqc_libcint_charges, only: ao_to_atom
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   implicit none

   real(dp), parameter :: A2B = 1.8897261254578281_dp
   integer, parameter :: NW = 3
   type(libcint_molecule_t) :: super, frag
   type(rhf_result_t) :: scf
   type(error_t) :: error
   type(direct_stats_t) :: stats
   real(dp), allocatable :: sb(:, :), fb(:, :), zero_s(:, :), zero_f(:, :)
   real(dp), allocatable :: d_tot(:, :), d_env(:, :), j_tot(:, :), j_env(:, :), j_own(:, :)
   real(dp), allocatable :: d_frag(:, :, :)
   integer, allocatable :: ao_atom(:), ao(:), owner(:), z(:)
   character(len=2), allocatable :: sym(:)
   real(dp), allocatable :: xyz(:, :)
   real(dp) :: mono(3, 3)
   real(dp) :: worst
   integer :: w, k, at, f, i
   character(len=160) :: line

   call logger%configure(info_level)
   mono = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                   0.0_dp, -0.7572_dp, 0.5865_dp, &
                   0.0_dp, 0.7572_dp, 0.5865_dp], [3, 3])
   allocate (z(3*NW), sym(3*NW), owner(3*NW), xyz(3, 3*NW))
   at = 0
   do w = 1, NW
      do k = 1, 3
         at = at + 1
         z(at) = merge(8, 1, k == 1)
         sym(at) = merge("O ", "H ", k == 1)
         xyz(:, at) = mono(:, k)
         xyz(3, at) = xyz(3, at) + real(w - 1, dp)*2.9_dp
         owner(at) = w
      end do
   end do
   xyz = xyz*A2B

   call build_libcint_molecule(z, sym, xyz, "6-31g", super, error)
   call schwarz_bounds(super, sb, error)
   call ao_to_atom(super, ao_atom)
   allocate (zero_s(super%nao, super%nao), source=0.0_dp)

   ! An isolated monomer density, replicated -- the identity does not care
   ! whether the densities are converged, only that they are densities.
   call build_libcint_molecule(z(1:3), sym(1:3), xyz(:, 1:3), "6-31g", frag, error)
   call schwarz_bounds(frag, fb, error)
   call run_libcint_rhf(frag, 10, 100, 1.0e-9_dp, 1.0e-7_dp, .false., scf, error)
   allocate (d_frag(frag%nao, frag%nao, NW))
   do f = 1, NW
      d_frag(:, :, f) = scf%density*(1.0_dp + 0.03_dp*real(f, dp))   !! deliberately unequal
   end do
   allocate (zero_f(frag%nao, frag%nao), source=0.0_dp)

   allocate (d_tot(super%nao, super%nao), source=0.0_dp)
   do f = 1, NW
      ao = pack([(i, i=1, super%nao)], owner(ao_atom) == f)
      d_tot(ao, ao) = d_frag(:, :, f)
   end do
   allocate (j_tot(super%nao, super%nao))
   call build_fock_direct(super, zero_s, d_tot, sb, j_tot, stats, error, &
                          k_scale=0.0_dp, j_scale=1.0_dp)

   worst = 0.0_dp
   allocate (j_env(super%nao, super%nao), j_own(frag%nao, frag%nao))
   do f = 1, NW
      ao = pack([(i, i=1, super%nao)], owner(ao_atom) == f)

      ! The slow, obvious way: J from the environment density alone.
      d_env = d_tot
      d_env(ao, ao) = 0.0_dp
      call build_fock_direct(super, zero_s, d_env, sb, j_env, stats, error, &
                             k_scale=0.0_dp, j_scale=1.0_dp)

      ! The way the module does it.
      call build_fock_direct(frag, zero_f, d_frag(:, :, f), fb, j_own, stats, error, &
                             k_scale=0.0_dp, j_scale=1.0_dp)

      write (line, "(a,i0,a,es12.4)") "   fragment ", f, "   largest difference ", &
         maxval(abs(j_env(ao, ao) - (j_tot(ao, ao) - j_own)))
      call logger%info(trim(line))
      worst = max(worst, maxval(abs(j_env(ao, ao) - (j_tot(ao, ao) - j_own))))
   end do

   write (line, "(a,es12.4)") "worst over all fragments: ", worst
   call logger%info(trim(line))
   if (worst > 1.0e-11_dp) then
      call logger%error("the subtraction is NOT the environment Coulomb operator")
      stop 1
   end if
   call logger%info("the fast embedding is the slow one")
end program probe_esp
