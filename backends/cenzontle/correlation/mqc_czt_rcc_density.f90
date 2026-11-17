!! The reduced densities of closed-shell CCSD in spatial orbitals
module mqc_czt_rcc_density
   !! What the amplitudes and the Lambda amplitudes together say about the
   !! electrons: the density a property is the trace of, and the density the
   !! gradient contracts with derivative integrals.
   !!
   !! These are the *correlation* parts. The reference determinant's own
   !! density -- two electrons in each occupied orbital -- is not included, so
   !! the one-particle density here has no trace and is what multiplies the
   !! Fock matrix in the correlation energy.
   !!
   !! **Equations.** PySCF's `cc/ccsd_rdm.py`, term for term, on arrays laid out
   !! as that module lays them out. The amplitude and Lambda equations next
   !! door follow the same package, so one set of conventions runs through all
   !! three.
   !!
   !! Unrelaxed: the orbitals are those of the reference and do not respond.
   !! The response is added by whoever needs it, because it depends on what is
   !! being differentiated.
   use pic_types, only: dp
   use mqc_tensor_contract, only: contract, transposed
   use mqc_czt_rcc, only: rcc_state_t
   use mqc_czt_rcc_lambda, only: rcc_lambda_t
   implicit none
   private

   public :: rccsd_one_particle_density

contains

   subroutine rccsd_one_particle_density(state, lambda, gamma)
      !! The correlation part of the spin-summed one-particle density
      !!
      !! Over the active orbitals, occupied first and then virtual, and
      !! symmetric. The occupied block traces to minus the virtual block's
      !! trace: correlation moves electrons and makes none.
      type(rcc_state_t), intent(in) :: state
      type(rcc_lambda_t), intent(in) :: lambda
      real(dp), allocatable, intent(out) :: gamma(:, :)   !! (no + nv, no + nv)

      real(dp), allocatable :: doo(:, :), dvv(:, :), dvo(:, :)
      real(dp), allocatable :: xtv(:, :), xt1(:, :), xt2(:, :), theta(:, :, :, :)
      integer :: no, nv, n

      no = state%no
      nv = state%nv
      n = no + nv
      associate (t1 => state%t1, t2 => state%t2, l1 => lambda%l1, l2 => lambda%l2)
         allocate (doo(no, no), dvv(nv, nv), xtv(no, no), xt1(no, no), xt2(nv, nv), &
                   source=0.0_dp)
         call contract("ja,ia->ij", t1, l1, doo, -1.0_dp)
         call contract("ia,ib->ab", t1, l1, dvv, 1.0_dp)
         call contract("ie,me->im", t1, l1, xtv, 1.0_dp)
         dvo = transpose(t1)
         call contract("im,ma->ai", xtv, t1, dvo, -1.0_dp)

         theta = 2.0_dp*t2 - transposed(t2, [1, 2, 4, 3])
         call contract("jkab,ikab->ij", theta, l2, doo, -1.0_dp)
         call contract("jica,jicb->ab", theta, l2, dvv, 1.0_dp)
         call contract("mnef,inef->mi", l2, theta, xt1, 1.0_dp)
         call contract("mnaf,mnef->ea", l2, theta, xt2, 1.0_dp)
         call contract("imae,me->ai", theta, l1, dvo, 1.0_dp)
         call contract("mi,ma->ai", xt1, t1, dvo, -1.0_dp)
         call contract("ie,ae->ai", t1, xt2, dvo, -1.0_dp)

         allocate (gamma(n, n))
         gamma(1:no, 1:no) = doo + transpose(doo)
         gamma(1:no, no + 1:n) = l1 + transpose(dvo)
         gamma(no + 1:n, 1:no) = transpose(gamma(1:no, no + 1:n))
         gamma(no + 1:n, no + 1:n) = dvv + transpose(dvv)
      end associate
   end subroutine rccsd_one_particle_density

end module mqc_czt_rcc_density
