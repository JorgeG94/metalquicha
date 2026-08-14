!! What an EFP interaction energy costs, and how it scales
!!
!!     OMP_NUM_THREADS=1  ./build/bench_efp
!!     OMP_NUM_THREADS=40 ./build/bench_efp
!!
!! **A dimer cannot be benchmarked.** Two waters are ten expansion points and 25
!! pairs; the whole energy is microseconds and any timing is noise plus the cost of
!! reading the potential. The loops only become the cost at cluster size, which is
!! also the only size anyone would want this fast for, so that is what this builds:
!! `n` copies of one water potential on a cubic lattice.
!!
!! Reported per term, because they scale differently. Electrostatics is a loop over
!! pairs of expansion points -- five per fragment -- so it grows as the square of
!! the fragment count. Polarization is a loop over pairs of polarizable points,
!! four per fragment, *inside* a self-consistent iteration, so it grows as the
!! square again but with an iteration count on top that itself depends on how
!! strongly the fragments couple.
!!
!! The same energy is printed alongside the timing, at every size and thread count.
!! A parallel reduction over floating point is not bit-reproducible -- the summation
!! order changes with the schedule -- so what is asserted is agreement to a
!! tolerance, not equality, and the tolerance is what tells you whether a
!! disagreement is rounding or a race.
program bench_efp
   use pic_types, only: dp, int64
   use mqc_error, only: error_t
   use mqc_efp_potential, only: efp_potential_t, make_efp_potential, &
                                write_efp_potential
   use mqc_efp_read, only: efp_fragment_t, read_efp_potential
   use mqc_efp_interaction, only: efp_system_t, build_efp_system, &
                                  electrostatic_energy, polarization_energy, &
                                  dispersion_energy_e6
   use omp_lib, only: omp_get_max_threads
   implicit none

   real(dp), parameter :: ANG = 1.0_dp/0.52917724924_dp
   !> Lattice spacing, in Angstrom. Loose enough that the induced dipoles converge
   !> and tight enough that nothing is negligible.
   real(dp), parameter :: SPACING = 3.5_dp
   integer, parameter :: SIZES(4) = [8, 27, 64, 125]

   type(efp_potential_t) :: pot
   type(efp_fragment_t) :: seed
   type(error_t) :: err
   real(dp) :: c(3, 3)
   integer :: z(3)
   integer :: k
   character(len=2) :: symbols(3)
   character(len=*), parameter :: path = "/tmp/bench_efp_water.efp"

   z = [8, 1, 1]
   symbols = ["O ", "H ", "H "]
   c = reshape([0.00000000000000_dp, 0.00000000009155_dp, 0.10077199490609_dp, &
                0.00000000000000_dp, 0.77250895271063_dp, -0.46780199741728_dp, &
                0.00000000000000_dp, -0.77250895280218_dp, -0.46780199748881_dp], &
               [3, 3])*ANG

   write (*, "(A,I0,A)") "  water lattice, ", omp_get_max_threads(), " threads"
   write (*, "(A)") ""

   call make_efp_potential(z, symbols, c, "6-31g*", "WATER", pot, err)
   if (err%has_error()) then
      write (*, "(A,A)") "  building the potential failed: ", err%get_message()
      error stop 1
   end if
   call write_efp_potential(pot, path, err)
   call read_efp_potential(path, seed, err)
   if (err%has_error()) then
      write (*, "(A,A)") "  reading it back failed: ", err%get_message()
      error stop 1
   end if

   write (*, "(A)") "  fragments   points     electrostatics        polarization"// &
      "         dispersion"
   do k = 1, size(SIZES)
      call one_size(SIZES(k), seed)
   end do
   write (*, "(A)") ""
   write (*, "(A)") "[bench] done"

   call seed%destroy()
   call pot%destroy()

contains

   subroutine one_size(n_frag, seed)
      integer, intent(in) :: n_frag
      type(efp_fragment_t), intent(in) :: seed

      type(efp_fragment_t), allocatable :: frags(:)
      type(efp_system_t) :: system
      type(error_t) :: err
      real(dp), allocatable :: translations(:, :)
      real(dp) :: t0, t1, e_stat, e_pol, e_disp
      real(dp) :: s_stat, s_pol, s_disp
      integer :: side, i, j, l, at

      side = nint(real(n_frag, dp)**(1.0_dp/3.0_dp))
      allocate (frags(n_frag), translations(3, n_frag))
      at = 0
      do i = 0, side - 1
         do j = 0, side - 1
            do l = 0, side - 1
               at = at + 1
               if (at > n_frag) exit
               frags(at) = seed
               translations(:, at) = [real(i, dp), real(j, dp), real(l, dp)]*SPACING*ANG
            end do
         end do
      end do

      call build_efp_system(frags, translations, system, err)
      if (err%has_error()) then
         write (*, "(A,A)") "  building the system failed: ", err%get_message()
         return
      end if

      call seconds(t0)
      e_stat = electrostatic_energy(system, 3, screen=.true.)
      call seconds(t1)
      s_stat = t1 - t0

      call seconds(t0)
      e_pol = polarization_energy(system, frags, err)
      call seconds(t1)
      s_pol = t1 - t0

      call seconds(t0)
      e_disp = dispersion_energy_e6(system, frags)
      call seconds(t1)
      s_disp = t1 - t0

      write (*, "(A,I7,I9,3(F12.6,A,F8.4,A))") "  ", n_frag, system%n_points, &
         e_stat, " (", s_stat, " s)", e_pol, " (", s_pol, " s)", &
         e_disp, " (", s_disp, " s)"

      call system%destroy()
      deallocate (frags, translations)
   end subroutine one_size

   subroutine seconds(t)
      !! Wall clock, which is what threading is supposed to reduce; `cpu_time` sums
      !! over threads and would report no speedup at all.
      real(dp), intent(out) :: t
      integer(int64) :: count, rate
      call system_clock(count, rate)
      t = real(count, dp)/real(rate, dp)
   end subroutine seconds

end program bench_efp
