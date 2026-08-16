program probe_esp
   !! Is a neighbour's contribution to the embedding really local?
   !!
   !! With a `resppc` cutoff, [[mqc_libcint_fmo]] builds the exact part of the
   !! embedding over a molecule holding only the fragment and its near
   !! neighbours -- not the whole system. That is what stops the cost per
   !! fragment growing once the system is bigger than a neighbourhood, and it is
   !! correct only if the Coulomb operator that neighbour K exerts on fragment X
   !! is unchanged by the existence of some third fragment L elsewhere.
   !!
   !! It should be: the integrals are (mn|ls) with mn on X and ls on K, and L's
   !! basis functions appear in neither. But "should be" is exactly the kind of
   !! assumption that is worth one direct measurement, because if it were wrong
   !! the failure would look like a slightly-off energy rather than anything
   !! obviously broken.
   !!
   !! So: three waters. Fragment 1's embedding from fragment 2 alone, computed
   !! once over the pair and once over all three, and the two compared.
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

   ! Fragment 1's AO block, and fragment 2's density alone in the supersystem.
   worst = 0.0_dp
   allocate (j_env(super%nao, super%nao), j_own(frag%nao, frag%nao))
   allocate (d_env(super%nao, super%nao))
   do f = 2, NW
      ao = pack([(i, i=1, super%nao)], owner(ao_atom) == 1)

      ! Over all three waters, with only fragment f carrying density.
      d_env = 0.0_dp
      block
         integer, allocatable :: ao_f(:)
         ao_f = pack([(i, i=1, super%nao)], owner(ao_atom) == f)
         d_env(ao_f, ao_f) = d_frag(:, :, f)
      end block
      call build_fock_direct(super, zero_s, d_env, sb, j_env, stats, error, &
                             k_scale=0.0_dp, j_scale=1.0_dp)

      ! Over the pair {1, f} alone, laid out the way the module lays it out:
      ! the embedded fragment first, its neighbour after.
      block
         type(libcint_molecule_t) :: local
         real(dp), allocatable :: lb(:, :), dl(:, :), zl(:, :), jl(:, :)
         integer, allocatable :: pair(:)
         integer :: n1, nf
         pair = [(i, i=1, 3), (i, i=3*(f - 1) + 1, 3*f)]
         call build_libcint_molecule(z(pair), sym(pair), xyz(:, pair), "6-31g", local, error)
         call schwarz_bounds(local, lb, error)
         n1 = frag%nao
         nf = frag%nao
         allocate (dl(local%nao, local%nao), source=0.0_dp)
         allocate (zl(local%nao, local%nao), source=0.0_dp)
         allocate (jl(local%nao, local%nao))
         dl(n1 + 1:n1 + nf, n1 + 1:n1 + nf) = d_frag(:, :, f)
         call build_fock_direct(local, zl, dl, lb, jl, stats, error, &
                                k_scale=0.0_dp, j_scale=1.0_dp)

         write (line, "(a,i0,a,es12.4)") "   fragment 1 embedded by fragment ", f, &
            "   largest difference ", maxval(abs(j_env(ao, ao) - jl(1:n1, 1:n1)))
         call logger%info(trim(line))
         worst = max(worst, maxval(abs(j_env(ao, ao) - jl(1:n1, 1:n1))))
      end block
   end do

   write (line, "(a,es12.4)") "worst over all fragments: ", worst
   call logger%info(trim(line))
   if (worst > 1.0e-11_dp) then
      call logger%error("a neighbour contribution is NOT local; the cutoff is unsound")
      stop 1
   end if
   call logger%info("a neighbour contribution is local, so the cutoff is sound")
end program probe_esp
