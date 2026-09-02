!! Periodic table data and element utilities
module mqc_elements
   !! Provides atomic numbers, element symbols, and atomic masses for the complete
   !! periodic table (elements 1-118) with conversion functions between representations.
   use pic_ascii, only: to_upper, to_lower
   use mqc_atomic_radii, only: covalent_radius_cordero
   use pic_types, only: dp
   implicit none
   private

   public :: element_symbol_to_number  !! Convert element symbol to atomic number
   public :: element_number_to_symbol  !! Convert atomic number to element symbol
   public :: element_mass              !! Get atomic mass by atomic number
   public :: element_covalent_radius   !! Get covalent radius by atomic number
   public :: element_vdw_radius        !! Get van der Waals radius by atomic number
   ! TODO: refactr to use findloc
   ! Periodic table data as module-level parameters
   integer, parameter :: n_elements = 118
   character(len=2), parameter :: element_symbols(n_elements) = [character(len=2) :: &
      !! Element symbols for the complete periodic table (H through Og)
      !! Ordered by atomic number from 1 to 118
                                                                 ! for some reason this is how the formatted formats this (????)
                                                                 "H", "He", &
                                                                 "Li", "Be", "B", "C", "N", "O", "F", "Ne", &
                                                                 "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", &
               "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", &
               "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", &
                   "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", &
                                "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", &
                    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", &
                                 "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]

   real(dp), parameter :: element_masses(n_elements) = [ &
      !! Standard atomic masses in atomic mass units (amu)
      !! Based on IUPAC standard atomic weights, ordered by atomic number
                          ! for some reason this is how the formatted formats this (????)
                          1.008_dp, 4.0026_dp, &                                                               ! H-He
                          6.94_dp, 9.0122_dp, 10.81_dp, 12.011_dp, 14.007_dp, 15.999_dp, 18.998_dp, 20.180_dp, &  ! Li-Ne
                          22.990_dp, 24.305_dp, 26.982_dp, 28.085_dp, 30.974_dp, 32.06_dp, 35.45_dp, 39.948_dp, &  ! Na-Ar
                          39.098_dp, 40.078_dp, 44.956_dp, 47.867_dp, 50.942_dp, 51.996_dp, 54.938_dp, 55.845_dp, &  ! K-Fe
                          58.933_dp, 58.693_dp, 63.546_dp, 65.38_dp, 69.723_dp, 72.630_dp, 74.922_dp, 78.971_dp, &  ! Co-Se
                          79.904_dp, 83.798_dp, &                                                              ! Br-Kr
                          85.468_dp, 87.62_dp, 88.906_dp, 91.224_dp, 92.906_dp, 95.95_dp, 98.0_dp, 101.07_dp, &  ! Rb-Ru
                          102.91_dp, 106.42_dp, 107.87_dp, 112.41_dp, 114.82_dp, 118.71_dp, 121.76_dp, 127.60_dp, &  ! Rh-Te
                          126.90_dp, 131.29_dp, &                                                              ! I-Xe
                          132.91_dp, 137.33_dp, 138.91_dp, 140.12_dp, 140.91_dp, 144.24_dp, 145.0_dp, 150.36_dp, &  ! Cs-Sm
                          151.96_dp, 157.25_dp, 158.93_dp, 162.50_dp, 164.93_dp, 167.26_dp, 168.93_dp, 173.05_dp, &  ! Eu-Yb
                          174.97_dp, 178.49_dp, 180.95_dp, 183.84_dp, 186.21_dp, 190.23_dp, 192.22_dp, 195.08_dp, &  ! Lu-Pt
                          196.97_dp, 200.59_dp, 204.38_dp, 207.2_dp, 208.98_dp, 209.0_dp, 210.0_dp, 222.0_dp, &  ! Au-Rn
                          223.0_dp, 226.0_dp, 227.0_dp, 232.04_dp, 231.04_dp, 238.03_dp, 237.0_dp, 244.0_dp, &  ! Fr-Pu
                          243.0_dp, 247.0_dp, 247.0_dp, 251.0_dp, 252.0_dp, 257.0_dp, 258.0_dp, 259.0_dp, &  ! Am-No
                          262.0_dp, 267.0_dp, 268.0_dp, 271.0_dp, 272.0_dp, 270.0_dp, 276.0_dp, 281.0_dp, &  ! Lr-Ds
                          280.0_dp, 285.0_dp, 284.0_dp, 289.0_dp, 288.0_dp, 293.0_dp, 294.0_dp, 294.0_dp]   ! Rg-Og

   !! Van der Waals radii in Angstrom, Bondi (1964) with Rowland and Taylor's
   !! hydrogen. Zero where Bondi tabulates none -- see `element_vdw_radius`.
   integer, parameter :: n_vdw = 96
   real(dp), parameter :: vdw_radii(n_vdw) = [ &
                          1.10_dp, 1.40_dp, 1.81_dp, 1.53_dp, 1.92_dp, 1.70_dp, 1.55_dp, 1.52_dp, &
                          1.47_dp, 1.54_dp, 2.27_dp, 1.73_dp, 1.84_dp, 2.10_dp, 1.80_dp, 1.80_dp, &
                          1.75_dp, 1.88_dp, 2.75_dp, 2.31_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, &
                          0.00_dp, 0.00_dp, 0.00_dp, 1.63_dp, 1.40_dp, 1.39_dp, 1.87_dp, 2.11_dp, &
                          1.85_dp, 1.90_dp, 1.85_dp, 2.02_dp, 3.03_dp, 2.49_dp, 0.00_dp, 0.00_dp, &
                          0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 1.63_dp, 1.72_dp, 1.58_dp, &
                          1.93_dp, 2.17_dp, 2.06_dp, 2.06_dp, 1.98_dp, 2.16_dp, 3.43_dp, 2.68_dp, &
                          0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, &
                          0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, &
                          0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 1.75_dp, 1.66_dp, 1.55_dp, &
                          1.96_dp, 2.02_dp, 2.07_dp, 1.97_dp, 2.02_dp, 2.20_dp, 3.48_dp, 2.83_dp, &
                          0.00_dp, 0.00_dp, 0.00_dp, 1.86_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp]

contains

   pure function element_symbol_to_number(symbol) result(atomic_number)
      !! Convert element symbol to atomic number
      !! Covers the complete periodic table (elements 1-118)
      character(len=*), intent(in) :: symbol
      integer :: atomic_number

      character(len=2) :: sym

      ! Normalize: uppercase first letter, lowercase second
      sym = adjustl(symbol)
      if (len_trim(sym) >= 1) sym(1:1) = to_upper(sym(1:1))
      if (len_trim(sym) >= 2) sym(2:2) = to_lower(sym(2:2))

      ! Search for symbol in table
      atomic_number = findloc(element_symbols, sym, dim=1)

   end function element_symbol_to_number

   pure function element_number_to_symbol(atomic_number) result(symbol)
      !! Convert atomic number to element symbol
      !! Covers the complete periodic table (elements 1-118)
      integer, intent(in) :: atomic_number
      character(len=2) :: symbol

      select case (atomic_number)
      case (1:118)
         symbol = element_symbols(atomic_number)
      case default
         symbol = "Xx"  ! Unknown
      end select

   end function element_number_to_symbol

   pure function element_mass(atomic_number) result(mass)
      !! Return atomic mass in atomic mass units (amu) for a given atomic number
      !! Uses standard atomic weights from IUPAC
      integer, intent(in) :: atomic_number
      real(dp) :: mass

      select case (atomic_number)
      case (1:118)
         mass = element_masses(atomic_number)
      case default
         mass = 0.0_dp  ! Unknown element
      end select

   end function element_mass

   pure function element_vdw_radius(atomic_number) result(radius)
      !! Van der Waals radius in Angstrom, or 0 where none is tabulated
      !!
      !! Bondi's set (J. Phys. Chem. 68, 441 (1964)) with Rowland and Taylor's
      !! revision for hydrogen, which is the usual modern choice.
      !!
      !! **These are a convention, and an electrostatic-potential charge fit
      !! depends on which one is used.** CHELPG excludes grid points inside these
      !! spheres, so a different radius set moves the points the fit sees and
      !! therefore the charges it returns. Comparing charges against another
      !! program means checking it uses these radii, not merely that it calls
      !! the method CHELPG.
      !!
      !! Zero for anything untabulated, for the same reason the covalent radii
      !! refuse: a caller has to decide what to do about an element nobody
      !! measured, and a made-up sphere would silently exclude or include points.
      integer, intent(in) :: atomic_number
      real(dp) :: radius

      select case (atomic_number)
      case (1:n_vdw)
         radius = vdw_radii(atomic_number)
      case default
         radius = 0.0_dp
      end select
   end function element_vdw_radius

   pure function element_covalent_radius(atomic_number) result(radius)
      !! Covalent radius in Angstrom, or 0 where none is tabulated
      !!
      !! The Cordero set, which is the one bond perception wants. Kept as a
      !! name on `mqc_elements` because that is where callers look for
      !! per-element data; the table itself lives in `mqc_atomic_radii`
      !! alongside the other parametrisations it must not be confused with.
      integer, intent(in) :: atomic_number
      real(dp) :: radius

      radius = covalent_radius_cordero(atomic_number)
   end function element_covalent_radius

end module mqc_elements
