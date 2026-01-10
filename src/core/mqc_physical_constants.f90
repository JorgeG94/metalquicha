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

end module mqc_physical_constants
