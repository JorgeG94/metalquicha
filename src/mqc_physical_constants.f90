!! Physical constants and unit conversion factors
module mqc_physical_constants
   !! Contains fundamental physical constants and unit conversion factors
   !! used throughout the metalquicha codebase.
   !!
   !! All values are in atomic units unless otherwise specified.
   !! Reference: CODATA 2018 recommended values where applicable.
   use pic_types, only: dp
   implicit none
   private

   !---------------------------------------------------------------------------
   ! Fundamental Constants
   !---------------------------------------------------------------------------

   !> Bohr radius in Angstrom (a0 = 0.529177... A)
   real(dp), parameter, public :: BOHR_TO_ANGSTROM = 0.52917721092_dp

   !> Angstrom to Bohr conversion
   real(dp), parameter, public :: ANGSTROM_TO_BOHR = 1.0_dp/BOHR_TO_ANGSTROM

   !> Atomic mass unit to atomic units of mass (electron masses)
   !> 1 amu = 1822.888 m_e
   real(dp), parameter, public :: AMU_TO_AU = 1822.888_dp

   !> Atomic units of mass to amu
   real(dp), parameter, public :: AU_TO_AMU = 1.0_dp/AMU_TO_AU

   !> Atomic mass unit to kg (CODATA 2018)
   real(dp), parameter, public :: AMU_TO_KG = 1.66053906660e-27_dp

   !---------------------------------------------------------------------------
   ! Vibrational Spectroscopy Conversions
   !---------------------------------------------------------------------------

   !> Conversion factor from atomic units (Hartree/Bohr^2/amu) to cm^-1
   !> Derived from: sqrt(Hartree/(Bohr^2 * amu)) -> s^-1 -> cm^-1
   real(dp), parameter, public :: AU_TO_CM1 = 2.642461e7_dp

   !> Conversion factor from atomic units (Hartree/Bohr^2) to mdyne/Angstrom
   !> 1 Hartree/Bohr^2 = 15.569141 mdyne/A
   real(dp), parameter, public :: AU_TO_MDYNE_ANG = 15.569141_dp

   !> Conversion factor from atomic units of dipole derivatives to km/mol (IR intensity)
   !> From IUPAC: A = (pi * N_A * |d_mu/dQ|^2) / (3 * 4 * pi * epsilon_0 * c^2)
   !> Reference: IUPAC, Quantities, Units and Symbols in Physical Chemistry (1993)
   real(dp), parameter, public :: AU_TO_KMMOL = 1.7770969e6_dp

   !---------------------------------------------------------------------------
   ! Dipole Moment Conversions
   !---------------------------------------------------------------------------

   !> Conversion from atomic units (e*Bohr) to Debye
   !> 1 e*a0 = 2.541746 Debye
   real(dp), parameter, public :: AU_TO_DEBYE = 2.541746_dp

   !> Conversion from Debye to atomic units
   real(dp), parameter, public :: DEBYE_TO_AU = 1.0_dp/AU_TO_DEBYE

   !---------------------------------------------------------------------------
   ! Energy Conversions
   !---------------------------------------------------------------------------

   !> Hartree to eV
   real(dp), parameter, public :: HARTREE_TO_EV = 27.211386245988_dp

   !> Hartree to kcal/mol
   real(dp), parameter, public :: HARTREE_TO_KCALMOL = 627.5094740631_dp

   !> Hartree to kJ/mol
   real(dp), parameter, public :: HARTREE_TO_KJMOL = 2625.4996394799_dp

   !> Hartree to cal/mol
   real(dp), parameter, public :: HARTREE_TO_CALMOL = 627.5094740631_dp*1000.0_dp

   !> Hartree to J/mol
   real(dp), parameter, public :: HARTREE_TO_JMOL = 2625.4996394799_dp*1000.0_dp

   !> Calorie to Joule (thermochemical calorie)
   real(dp), parameter, public :: CAL_TO_J = 4.184_dp

   !---------------------------------------------------------------------------
   ! Thermochemistry Constants (CODATA 2018)
   !---------------------------------------------------------------------------

   !> Boltzmann constant in Hartree/K
   !> k_B = 1.380649e-23 J/K, 1 Hartree = 4.3597447222071e-18 J
   real(dp), parameter, public :: KB_HARTREE = 3.1668115634556e-6_dp

   !> Boltzmann constant in J/K (CODATA 2018, exact)
   real(dp), parameter, public :: KB_SI = 1.380649e-23_dp

   !> Planck constant in Hartree*s
   !> h = 6.62607015e-34 J*s
   real(dp), parameter, public :: H_HARTREE_S = 1.5198298460574e-16_dp

   !> Planck constant in J*s (CODATA 2018, exact)
   real(dp), parameter, public :: H_SI = 6.62607015e-34_dp

   !> Speed of light in cm/s
   real(dp), parameter, public :: C_CM_S = 2.99792458e10_dp

   !> cm^-1 to Kelvin conversion factor: theta_vib = (h*c/k_B) * nu
   !> This is h*c/k_B in cm (multiply by frequency in cm^-1 to get K)
   real(dp), parameter, public :: CM1_TO_KELVIN = 1.4387773538277_dp

   !> Gas constant R in cal/(mol*K) for thermochemistry output
   !> R = 1.98720425864 cal/(mol*K)
   real(dp), parameter, public :: R_CALMOLK = 1.98720425864_dp

   !> Gas constant R in Hartree/(mol*K)
   !> R = N_A * k_B
   real(dp), parameter, public :: R_HARTREE = 3.1668115634556e-6_dp

   !> Pressure: 1 atm in atomic units (Hartree/Bohr^3)
   !> 1 atm = 101325 Pa, 1 Bohr = 5.29177e-11 m, 1 Hartree = 4.3597e-18 J
   real(dp), parameter, public :: ATM_TO_AU = 3.39893097e-9_dp

   !> Pressure: 1 atm in Pa (exact by definition)
   real(dp), parameter, public :: ATM_TO_PA = 101325.0_dp

   !> Pi constant
   real(dp), parameter, public :: PI = 3.14159265358979323846_dp

   !> Avogadro's number (for reference, not directly used in atomic unit calculations)
   real(dp), parameter, public :: AVOGADRO = 6.02214076e23_dp

   !> Rotational constant conversion: amu*Angstrom^2 to GHz
   !> B = h / (8*pi^2*I) where I is in SI units
   !> For I in amu*Angstrom^2: B(GHz) = 505379.07 / I
   real(dp), parameter, public :: ROTCONST_AMUA2_TO_GHZ = 505379.07_dp

   !> Rotational temperature conversion: amu*Angstrom^2 to Kelvin
   !> theta_rot = h^2 / (8*pi^2*I*k_B)
   !> For I in amu*Angstrom^2: theta_rot(K) = 24.2637 / I
   real(dp), parameter, public :: ROTTEMP_AMUA2_TO_K = 24.2637_dp

   !> Classical limit for vibrational modes (u = theta_v/T)
   !> When u > this value, vibrational modes are considered frozen out
   real(dp), parameter, public :: VIB_CLASSICAL_LIMIT = 100.0_dp

end module mqc_physical_constants
