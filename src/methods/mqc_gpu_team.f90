!! The set of ranks that share one system's GPU work
module mqc_gpu_team
   !! Process-wide state for the multi-GPU mode: several MPI ranks, each bound
   !! to its own GPU, cooperating on a *single* system rather than on separate
   !! fragments.
   !!
   !! This is deliberately not the fragmented parallelism. There, one rank owns
   !! a whole subsystem and nothing has to be reduced mid-SCF; here every rank
   !! runs the same SCF and owns a slice of one quantity, so the ranks are in
   !! lockstep and a reduction happens inside every iteration.
   !!
   !! What is sliced is the exchange-correlation quadrature, and only that.
   !! cuEST 0.2 exposes no way to split a density-fitted J or K build -- the
   !! pair list is built from a threshold rather than handed a subset, and
   !! splitting the auxiliary basis would change the fitting metric rather than
   !! partition it -- so every rank still builds a full DF plan. The mode
   !! divides XC time and XC memory by the team size and nothing else. See
   !! `backends/cuest/CUEST.md` for why, and for the knobs that do reduce the
   !! DF plan.
   !!
   !! The team lives here, in `src/`, and knows nothing about cuEST or CUDA,
   !! for the same reason `mqc_cuest_iface` does: the backend cannot be part of
   !! the fpm build, and the driver that enables the team cannot reach it.
   !! Reductions are over host arrays; staging a device buffer through the host
   !! is the backend's business, because only the backend knows what a device
   !! pointer is.
   use pic_types, only: dp
   use pic_mpi_lib, only: comm_t, allreduce
   implicit none
   private

   public :: gpu_team_enable      !! Put this rank into a team
   public :: gpu_team_disable     !! Leave the team; reductions become no-ops
   public :: gpu_team_active      !! Whether a team of more than one rank is live
   public :: gpu_team_rank        !! This rank's index in the team, 0-based
   public :: gpu_team_size        !! How many ranks share the system
   public :: gpu_team_reduce      !! Sum a host array across the team, in place
   public :: gpu_team_reduce_scalar  !! Sum a scalar across the team, in place
   public :: gpu_team_agree       !! Make a logical unanimous across the team

   logical, save :: team_live = .false.
   integer, save :: team_rank = 0
   integer, save :: team_size = 1
   type(comm_t), save :: team_comm

contains

   subroutine gpu_team_enable(comm)
      !! Enrol this rank in the team spanning `comm`
      !!
      !! A single-rank communicator enables nothing: the mode would then cost a
      !! reduction per iteration and buy no division of work, and `active`
      !! staying false keeps every reduction below a no-op rather than a
      !! zero-length collective.
      type(comm_t), intent(in) :: comm

      team_comm = comm
      team_rank = comm%rank()
      team_size = comm%size()
      team_live = (team_size > 1)
   end subroutine gpu_team_enable

   subroutine gpu_team_disable()
      !! Leave the team
      team_live = .false.
      team_rank = 0
      team_size = 1
   end subroutine gpu_team_disable

   pure function gpu_team_active() result(active)
      !! Whether work should be sliced and results reduced
      logical :: active

      active = team_live
   end function gpu_team_active

   pure function gpu_team_rank() result(rank)
      !! This rank's slice index, 0-based; 0 when no team is live
      integer :: rank

      rank = team_rank
   end function gpu_team_rank

   pure function gpu_team_size() result(n_slices)
      !! Number of slices; 1 when no team is live
      integer :: n_slices

      n_slices = team_size
   end function gpu_team_size

   subroutine gpu_team_reduce(buffer)
      !! Sum a host array over the team, in place
      !!
      !! A no-op outside a team, so callers need no conditional of their own.
      real(dp), intent(inout) :: buffer(:)

      if (.not. team_live) return
      if (size(buffer) == 0) return
      call allreduce(team_comm, buffer)
   end subroutine gpu_team_reduce

   subroutine gpu_team_reduce_scalar(value)
      !! Sum a scalar over the team, in place
      real(dp), intent(inout) :: value

      if (.not. team_live) return
      call allreduce(team_comm, value)
   end subroutine gpu_team_reduce_scalar

   subroutine gpu_team_agree(flag)
      !! Force every rank to the same answer, true only if all ranks say true
      !!
      !! Used on the SCF's convergence test. The ranks should already agree --
      !! every term but the reduced XC one is computed redundantly from
      !! identical inputs, and cuEST is pinned to FP64 -- but "should" is doing
      !! a lot of work there, and the failure mode if it is wrong is not a
      !! wrong number, it is one rank leaving the loop while the others wait
      !! forever in the next reduction. One collective per iteration is cheap
      !! insurance against a hang.
      logical, intent(inout) :: flag

      real(dp) :: vote

      if (.not. team_live) return
      vote = 0.0_dp
      if (flag) vote = 1.0_dp
      call allreduce(team_comm, vote)
      flag = (vote >= real(team_size, dp) - 0.5_dp)
   end subroutine gpu_team_agree

end module mqc_gpu_team
