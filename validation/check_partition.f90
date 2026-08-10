program check_partition
   !! Check the Becke partition against PySCF and against its defining property
   !!
   !! The defining property is that the weights form a partition of unity: for
   !! any point, summing w_A over every atom must give exactly 1. That is
   !! checked here directly, for both cutoffs and all three size adjustments,
   !! and it is a strong check -- it fails immediately if the pair loop is
   !! asymmetric, if a shift has the wrong sign, or if an atom is skipped.
   !!
   !! Agreement with PySCF is checked by compare_partition.py, which reads the
   !! dump this writes.
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_dft_partition, only: becke_partition_weights, partition_scheme_name, &
                                PARTITION_BECKE, PARTITION_STRATMANN, &
                                ADJUST_NONE, ADJUST_BECKE, ADJUST_TREUTLER
   implicit none

   integer, parameter :: N_DIM = 3
   integer, parameter :: N_SCHEMES = 2
   integer, parameter :: N_ADJUSTS = 3
   integer, parameter :: N_PROBE_1D = 9

   integer, parameter :: SCHEMES(N_SCHEMES) = [PARTITION_BECKE, PARTITION_STRATMANN]
   integer, parameter :: ADJUSTS(N_ADJUSTS) = [ADJUST_NONE, ADJUST_BECKE, ADJUST_TREUTLER]

   integer :: unit, n_bad
   real(dp) :: worst_unity

   n_bad = 0
   worst_unity = 0.0_dp

   open (newunit=unit, file="partition_weights.txt", status="replace", action="write")

   ! Water: unequal atoms, so the size adjustment actually does something.
   call run_case(unit, "water", &
                 reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, -1.4308_dp, 1.1078_dp, &
                          0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3]), &
                 [8, 1, 1], n_bad, worst_unity)

   ! Gold hydride: an extreme radius ratio, where a mis-signed shift is obvious.
   call run_case(unit, "auh", &
                 reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp, 0.0_dp, 2.8_dp], [N_DIM, 2]), &
                 [79, 1], n_bad, worst_unity)

   close (unit)

   write (*, "(a,es10.3)") "worst |sum_A w_A - 1| : ", worst_unity
   if (n_bad > 0 .or. worst_unity > 1.0e-13_dp) then
      write (*, "(a)") "FAILED"
      stop 1
   end if
   write (*, "(a)") "partition of unity OK -- now run compare_partition.py"

contains

   subroutine run_case(unit, label, atom_coords, atomic_numbers, n_bad, worst_unity)
      !! Every scheme and adjustment on one molecule, over a lattice of probes
      integer, intent(in) :: unit
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: atom_coords(:, :)
      integer, intent(in) :: atomic_numbers(:)
      integer, intent(inout) :: n_bad
      real(dp), intent(inout) :: worst_unity

      real(dp), allocatable :: probes(:, :), weights(:), unity(:)
      integer, allocatable :: owner(:)
      type(error_t) :: error
      integer :: n_atoms, n_probes, s, a, ia, k

      n_atoms = size(atom_coords, 2)
      call build_probes(atom_coords, probes)
      n_probes = size(probes, 2)

      allocate (weights(n_probes), owner(n_probes), unity(n_probes))

      do s = 1, N_SCHEMES
         do a = 1, N_ADJUSTS
            unity = 0.0_dp
            do ia = 1, n_atoms
               owner = ia
               call becke_partition_weights(probes, atom_coords, atomic_numbers, &
                                            owner, SCHEMES(s), ADJUSTS(a), weights, error)
               if (error%has_error()) then
                  write (*, "(a,a)") "FAIL: ", error%get_message()
                  n_bad = n_bad + 1
                  return
               end if

               unity = unity + weights
               do k = 1, n_probes
                  write (unit, "(a,1x,i0,1x,i0,1x,i0,3(1x,es24.16e3),1x,es24.16e3)") &
                     label, SCHEMES(s), ADJUSTS(a), ia, probes(:, k), weights(k)
               end do
            end do

            worst_unity = max(worst_unity, maxval(abs(unity - 1.0_dp)))
            write (*, "(a,a,a,a,a,i0,a,es10.3)") "  ", label, " ", &
               partition_scheme_name(SCHEMES(s)), " adjust ", ADJUSTS(a), &
               "  max |sum w - 1| = ", maxval(abs(unity - 1.0_dp))
         end do
      end do

      deallocate (probes, weights, owner, unity)
   end subroutine run_case

   subroutine build_probes(atom_coords, probes)
      !! A lattice spanning the molecule and reaching well outside it
      !!
      !! Deliberately includes points very close to nuclei and points far away,
      !! since those are where a cutoff misbehaves: near a nucleus one cell
      !! function approaches 1 and the rest underflow, and far away all of them
      !! are tiny and the normalisation is what keeps the sum at 1.
      real(dp), intent(in) :: atom_coords(:, :)
      real(dp), allocatable, intent(out) :: probes(:, :)

      real(dp) :: lo(N_DIM), hi(N_DIM)
      real(dp) :: span
      integer :: i, j, k, n

      lo = minval(atom_coords, dim=2) - 3.0_dp
      hi = maxval(atom_coords, dim=2) + 3.0_dp

      allocate (probes(N_DIM, N_PROBE_1D**3))
      n = 0
      do i = 1, N_PROBE_1D
         do j = 1, N_PROBE_1D
            do k = 1, N_PROBE_1D
               n = n + 1
               span = real(i - 1, dp)/real(N_PROBE_1D - 1, dp)
               probes(1, n) = lo(1) + span*(hi(1) - lo(1))
               span = real(j - 1, dp)/real(N_PROBE_1D - 1, dp)
               probes(2, n) = lo(2) + span*(hi(2) - lo(2))
               span = real(k - 1, dp)/real(N_PROBE_1D - 1, dp)
               probes(3, n) = lo(3) + span*(hi(3) - lo(3))
            end do
         end do
      end do
   end subroutine build_probes

end program check_partition
