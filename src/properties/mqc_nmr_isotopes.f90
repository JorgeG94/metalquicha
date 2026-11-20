!! The isotope an NMR experiment observes for each element, and its g factor
module mqc_nmr_isotopes
   !! A spin-spin coupling in hertz belongs to a pair of isotopes, not of
   !! elements: it is the reduced coupling times the two nuclear g factors.
   !! This is the table that turns one into the other, for the isotope a
   !! spectrometer is usually tuned to. Carbon is 13C, not the abundant and
   !! spinless 12C.
   use pic_types, only: dp
   implicit none
   private

   public :: nmr_isotope
   public :: reduced_to_hertz
   public :: reduced_to_si

   integer, parameter :: N_TABLE = 17
   integer, parameter :: LABEL_LEN = 5
   character(len=LABEL_LEN), parameter :: LABELS(N_TABLE) = [character(len=LABEL_LEN) :: &
      "1H", "3He", "7Li", "9Be", "11B", "13C", "15N", "17O", "19F", "21Ne", "23Na", "25Mg", &
      "27Al", "29Si", "31P", "33S", "35Cl"]
   real(dp), parameter :: G_FACTORS(N_TABLE) = [ &
                          5.5856946893_dp, -4.255250615_dp, 2.170951_dp, -0.784993_dp, &
                          1.792433_dp, 1.4048236_dp, -0.56637768_dp, -0.757516_dp, &
                          5.257736_dp, -0.441198_dp, 1.478348_dp, -0.34218_dp, &
                          1.4566028_dp, -1.11058_dp, 2.26320_dp, 0.429214_dp, 0.5479162_dp]
      !! Magnetic moment in nuclear magnetons over the spin, H through Cl

   real(dp), parameter :: HARTREE_TO_HZ = 6.579683920502e15_dp
   real(dp), parameter :: HARTREE_TO_JOULE = 4.3597447222071e-18_dp
   real(dp), parameter :: NUCLEAR_MAGNETON_SI = 5.0507837461e-27_dp   ! J/T
   real(dp), parameter :: PROTON_ELECTRON_MASS = 1836.15267343_dp
   real(dp), parameter :: NUCLEAR_MAGNETON_AU = 0.5_dp/PROTON_ELECTRON_MASS

contains

   pure subroutine nmr_isotope(atomic_number, label, g_factor)
      !! The observed isotope of an element and its g factor
      !!
      !! `label` is blank and `g_factor` zero for an element beyond the
      !! table, whose couplings are then reported reduced only.
      integer, intent(in) :: atomic_number
      character(len=*), intent(out) :: label
      real(dp), intent(out) :: g_factor

      label = ""
      g_factor = 0.0_dp
      if (atomic_number < 1 .or. atomic_number > N_TABLE) return
      label = trim(LABELS(atomic_number))
      g_factor = G_FACTORS(atomic_number)
   end subroutine nmr_isotope

   pure function reduced_to_hertz(reduced, g_a, g_b) result(j)
      !! `J` in Hz from a reduced coupling in atomic units and two g factors
      real(dp), intent(in) :: reduced, g_a, g_b
      real(dp) :: j

      j = HARTREE_TO_HZ*NUCLEAR_MAGNETON_AU**2*g_a*g_b*reduced
   end function reduced_to_hertz

   pure function reduced_to_si(reduced) result(k)
      !! The reduced coupling `K` in 10^19 T^2/J, which no isotope enters
      real(dp), intent(in) :: reduced
      real(dp) :: k

      k = HARTREE_TO_JOULE*NUCLEAR_MAGNETON_AU**2*reduced/NUCLEAR_MAGNETON_SI**2*1.0e-19_dp
   end function reduced_to_si

end module mqc_nmr_isotopes
