!! Thermochemistry calculations using the Rigid Rotor Harmonic Oscillator (RRHO) approximation
module mqc_thermochemistry
   !! Computes thermodynamic properties from vibrational frequencies and molecular geometry.
   !!
   !! This module implements standard RRHO thermochemistry calculations including:
   !! - Zero-point vibrational energy (ZPE)
   !! - Translational, rotational, vibrational, and electronic contributions
   !! - Thermal corrections to energy, enthalpy, and Gibbs free energy
   !!
   !! Default conditions: T = 298.15 K, P = 1 atm
   !! Output follows Gaussian-style formatting for compatibility.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_physical_constants, only: BOHR_TO_ANGSTROM, AMU_TO_AU, AMU_TO_KG, &
                                     KB_HARTREE, KB_SI, H_SI, CM1_TO_KELVIN, &
                                     R_CALMOLK, R_HARTREE, ATM_TO_PA, CAL_TO_J, &
                                     PI, ROTCONST_AMUA2_TO_GHZ, ROTTEMP_AMUA2_TO_K, &
                                     VIB_CLASSICAL_LIMIT, &
                                     HARTREE_TO_KCALMOL, HARTREE_TO_JMOL, HARTREE_TO_CALMOL
   use mqc_elements, only: element_mass
   use pic_lapack_interfaces, only: pic_syev
   implicit none
   private

   public :: thermochemistry_result_t
   public :: compute_moments_of_inertia
   public :: compute_rotational_constants
   public :: compute_zpe
   public :: compute_translational_thermo
   public :: compute_rotational_thermo
   public :: compute_vibrational_thermo
   public :: compute_electronic_entropy
   public :: compute_thermochemistry
   public :: print_thermochemistry

   !> Default temperature in Kelvin
   real(dp), parameter, public :: DEFAULT_TEMPERATURE = 298.15_dp

   !> Default pressure in atm
   real(dp), parameter, public :: DEFAULT_PRESSURE = 1.0_dp

   !> Default symmetry number
   integer, parameter, public :: DEFAULT_SYMMETRY_NUMBER = 1

   !> Default spin multiplicity (singlet)
   integer, parameter, public :: DEFAULT_SPIN_MULTIPLICITY = 1

   !> Threshold for considering a moment of inertia as zero (linear molecule detection)
   real(dp), parameter :: LINEAR_THRESHOLD = 1.0e-6_dp

   !> Threshold for considering a frequency as imaginary
   real(dp), parameter :: IMAG_FREQ_THRESHOLD = 0.0_dp

   type :: thermochemistry_result_t
      !! Container for thermochemistry calculation results
      real(dp) :: temperature = DEFAULT_TEMPERATURE     !! Temperature in K
      real(dp) :: pressure = DEFAULT_PRESSURE           !! Pressure in atm
      integer :: symmetry_number = DEFAULT_SYMMETRY_NUMBER  !! Rotational symmetry number
      integer :: spin_multiplicity = DEFAULT_SPIN_MULTIPLICITY  !! Electronic spin multiplicity
      logical :: is_linear = .false.                    !! True if molecule is linear

      ! Molecular properties
      real(dp) :: total_mass = 0.0_dp                   !! Total mass in amu
      real(dp) :: moments(3) = 0.0_dp                   !! Principal moments of inertia in amu*Angstrom^2
      real(dp) :: rot_const(3) = 0.0_dp                 !! Rotational constants in GHz

      ! Zero-point energy
      real(dp) :: zpe_hartree = 0.0_dp                  !! ZPE in Hartree
      real(dp) :: zpe_kcalmol = 0.0_dp                  !! ZPE in kcal/mol

      ! Thermal energy contributions (Hartree)
      real(dp) :: E_trans = 0.0_dp                      !! Translational thermal energy
      real(dp) :: E_rot = 0.0_dp                        !! Rotational thermal energy
      real(dp) :: E_vib = 0.0_dp                        !! Vibrational thermal energy (excluding ZPE)
      real(dp) :: E_elec = 0.0_dp                       !! Electronic thermal energy (always 0)

      ! Entropy contributions (cal/(mol*K))
      real(dp) :: S_trans = 0.0_dp                      !! Translational entropy
      real(dp) :: S_rot = 0.0_dp                        !! Rotational entropy
      real(dp) :: S_vib = 0.0_dp                        !! Vibrational entropy
      real(dp) :: S_elec = 0.0_dp                       !! Electronic entropy

      ! Heat capacity contributions (cal/(mol*K))
      real(dp) :: Cv_trans = 0.0_dp                     !! Translational heat capacity
      real(dp) :: Cv_rot = 0.0_dp                       !! Rotational heat capacity
      real(dp) :: Cv_vib = 0.0_dp                       !! Vibrational heat capacity

      ! Summary quantities (Hartree)
      real(dp) :: thermal_correction_energy = 0.0_dp    !! E_tot = ZPE + E_trans + E_rot + E_vib
      real(dp) :: thermal_correction_enthalpy = 0.0_dp  !! H = E_tot + RT
      real(dp) :: thermal_correction_gibbs = 0.0_dp     !! G = H - TS

      ! Partition functions
      real(dp) :: q_trans = 0.0_dp                      !! Translational partition function
      real(dp) :: q_rot = 0.0_dp                        !! Rotational partition function
      real(dp) :: q_vib = 1.0_dp                        !! Vibrational partition function

      ! Counts
      integer :: n_real_freqs = 0                       !! Number of real vibrational frequencies
      integer :: n_imag_freqs = 0                       !! Number of imaginary frequencies (skipped)
   end type thermochemistry_result_t

contains

   subroutine compute_moments_of_inertia(coords, atomic_numbers, n_atoms, &
                                         center_of_mass, moments, principal_axes, is_linear, total_mass)
      !! Compute the principal moments of inertia and detect linear molecules.
      !!
      !! Calculates the center of mass, inertia tensor, and diagonalizes to get
      !! principal moments. A molecule is considered linear if one moment is
      !! essentially zero (< LINEAR_THRESHOLD).
      real(dp), intent(in) :: coords(:, :)       !! Atomic coordinates (3, n_atoms) in Bohr
      integer, intent(in) :: atomic_numbers(:)   !! Atomic numbers
      integer, intent(in) :: n_atoms             !! Number of atoms
      real(dp), intent(out) :: center_of_mass(3)  !! Center of mass in Angstrom
      real(dp), intent(out) :: moments(3)        !! Principal moments in amu*Angstrom^2
      real(dp), intent(out) :: principal_axes(3, 3)  !! Principal axis vectors (columns)
      logical, intent(out) :: is_linear          !! True if molecule is linear
      real(dp), intent(out) :: total_mass        !! Total mass in amu

      real(dp) :: coords_ang(3, n_atoms)
      real(dp) :: masses(n_atoms)
      real(dp) :: rel_coords(3, n_atoms)
      real(dp) :: inertia_tensor(3, 3)
      integer :: i, info
      real(dp) :: mass_i, x, y, z

      ! Convert coordinates to Angstrom
      coords_ang = coords*BOHR_TO_ANGSTROM

      ! Get atomic masses
      total_mass = 0.0_dp
      do i = 1, n_atoms
         masses(i) = element_mass(atomic_numbers(i))
         total_mass = total_mass + masses(i)
      end do

      ! Compute center of mass
      center_of_mass = 0.0_dp
      do i = 1, n_atoms
         center_of_mass(:) = center_of_mass(:) + masses(i)*coords_ang(:, i)
      end do
      center_of_mass = center_of_mass/total_mass

      ! Compute coordinates relative to center of mass
      do i = 1, n_atoms
         rel_coords(:, i) = coords_ang(:, i) - center_of_mass(:)
      end do

      ! Build inertia tensor
      inertia_tensor = 0.0_dp
      do i = 1, n_atoms
         mass_i = masses(i)
         x = rel_coords(1, i)
         y = rel_coords(2, i)
         z = rel_coords(3, i)

         ! Diagonal elements: I_xx = sum(m * (y^2 + z^2)), etc.
         inertia_tensor(1, 1) = inertia_tensor(1, 1) + mass_i*(y*y + z*z)
         inertia_tensor(2, 2) = inertia_tensor(2, 2) + mass_i*(x*x + z*z)
         inertia_tensor(3, 3) = inertia_tensor(3, 3) + mass_i*(x*x + y*y)

         ! Off-diagonal elements: I_xy = -sum(m * x * y), etc.
         inertia_tensor(1, 2) = inertia_tensor(1, 2) - mass_i*x*y
         inertia_tensor(1, 3) = inertia_tensor(1, 3) - mass_i*x*z
         inertia_tensor(2, 3) = inertia_tensor(2, 3) - mass_i*y*z
      end do

      ! Symmetrize
      inertia_tensor(2, 1) = inertia_tensor(1, 2)
      inertia_tensor(3, 1) = inertia_tensor(1, 3)
      inertia_tensor(3, 2) = inertia_tensor(2, 3)

      ! Diagonalize to get principal moments
      principal_axes = inertia_tensor
      call pic_syev(principal_axes, moments, 'V', 'U', info)

      if (info /= 0) then
         call logger%warning("Failed to diagonalize inertia tensor, info = "// &
                             trim(adjustl(to_char(info))))
         moments = 0.0_dp
         is_linear = .false.
         return
      end if

      ! Check for linear molecule: one moment should be ~0
      ! Moments are returned in ascending order
      is_linear = (moments(1) < LINEAR_THRESHOLD)

   end subroutine compute_moments_of_inertia

   subroutine compute_rotational_constants(moments, is_linear, rot_const)
      !! Convert moments of inertia to rotational constants in GHz.
      !!
      !! B = h / (8 * pi^2 * I) where I is in SI units.
      !! For I in amu*Angstrom^2: B(GHz) = 505379.07 / I
      real(dp), intent(in) :: moments(3)      !! Moments in amu*Angstrom^2
      logical, intent(in) :: is_linear        !! True if linear molecule
      real(dp), intent(out) :: rot_const(3)   !! Rotational constants in GHz

      integer :: i

      rot_const = 0.0_dp

      if (is_linear) then
         ! For linear molecules, only one rotational constant matters
         ! Use the largest moment (moments are sorted ascending)
         if (moments(3) > LINEAR_THRESHOLD) then
            rot_const(1) = ROTCONST_AMUA2_TO_GHZ/moments(3)
         end if
      else
         ! For nonlinear molecules, compute all three
         do i = 1, 3
            if (moments(i) > LINEAR_THRESHOLD) then
               rot_const(i) = ROTCONST_AMUA2_TO_GHZ/moments(i)
            end if
         end do
      end if

   end subroutine compute_rotational_constants

   subroutine compute_zpe(frequencies, n_freqs, n_real, zpe_hartree, zpe_kcalmol)
      !! Compute zero-point vibrational energy from frequencies.
      !!
      !! ZPE = (1/2) * h * sum(nu_i) for all real frequencies.
      !! Imaginary frequencies (negative values) are skipped with a warning.
      real(dp), intent(in) :: frequencies(:)  !! Vibrational frequencies in cm^-1
      integer, intent(in) :: n_freqs          !! Total number of frequencies
      integer, intent(out) :: n_real          !! Number of real (positive) frequencies used
      real(dp), intent(out) :: zpe_hartree    !! ZPE in Hartree
      real(dp), intent(out) :: zpe_kcalmol    !! ZPE in kcal/mol

      integer :: i
      real(dp) :: freq_sum
      integer :: n_imag

      freq_sum = 0.0_dp
      n_real = 0
      n_imag = 0

      do i = 1, n_freqs
         if (frequencies(i) > IMAG_FREQ_THRESHOLD) then
            freq_sum = freq_sum + frequencies(i)
            n_real = n_real + 1
         else if (frequencies(i) < IMAG_FREQ_THRESHOLD) then
            n_imag = n_imag + 1
         end if
         ! Frequencies exactly at threshold (typically trans/rot modes ~0) are skipped
      end do

      if (n_imag > 0) then
         call logger%warning("Thermochemistry: "//trim(adjustl(to_char(n_imag)))// &
                             " imaginary frequency(ies) skipped")
      end if

      ! ZPE = 0.5 * sum(h * c * nu) where nu is in cm^-1
      ! In atomic units: ZPE = 0.5 * sum(nu * CM1_TO_KELVIN * KB_HARTREE)
      ! Actually: h*c*nu [cm^-1] = h*c*nu [J] = nu * CM1_TO_KELVIN * k_B [J]
      ! So ZPE [Hartree] = 0.5 * sum(nu) * CM1_TO_KELVIN * KB_HARTREE
      zpe_hartree = 0.5_dp*freq_sum*CM1_TO_KELVIN*KB_HARTREE
      zpe_kcalmol = zpe_hartree*HARTREE_TO_KCALMOL

   end subroutine compute_zpe

   subroutine compute_translational_thermo(total_mass, temperature, pressure, E, S, Cv)
      !! Compute translational contributions to thermodynamic properties.
      !!
      !! Uses ideal gas partition function (Sackur-Tetrode equation for entropy).
      !! E_trans = 3/2 * R * T
      !! S_trans = R * [5/2 + ln((2*pi*m*k*T/h^2)^(3/2) * k*T/P)]
      !! Cv_trans = 3/2 * R
      real(dp), intent(in) :: total_mass   !! Total mass in amu
      real(dp), intent(in) :: temperature  !! Temperature in K
      real(dp), intent(in) :: pressure     !! Pressure in atm
      real(dp), intent(out) :: E           !! Thermal energy in Hartree
      real(dp), intent(out) :: S           !! Entropy in cal/(mol*K)
      real(dp), intent(out) :: Cv          !! Heat capacity in cal/(mol*K)

      real(dp) :: mass_kg, T, P_pa
      real(dp) :: lambda_cubed, V_molar, qt

      ! Convert inputs to SI
      mass_kg = total_mass*AMU_TO_KG
      T = temperature
      P_pa = pressure*ATM_TO_PA

      ! Thermal de Broglie wavelength cubed: lambda^3 = (h^2 / (2*pi*m*k*T))^(3/2)
      lambda_cubed = (H_SI*H_SI/(2.0_dp*PI*mass_kg*KB_SI*T))**1.5_dp

      ! Molar volume at given T and P (ideal gas): V = R*T/P = k*T/P per molecule
      V_molar = KB_SI*T/P_pa  ! m^3 per molecule

      ! Translational partition function per molecule
      qt = V_molar/lambda_cubed

      ! Thermal energy: E = 3/2 * R * T
      E = 1.5_dp*R_HARTREE*temperature

      ! Entropy (Sackur-Tetrode): S = R * [5/2 + ln(qt)]
      S = R_CALMOLK*(2.5_dp + log(qt))

      ! Heat capacity: Cv = 3/2 * R
      Cv = 1.5_dp*R_CALMOLK

   end subroutine compute_translational_thermo

   subroutine compute_rotational_thermo(moments, temperature, symmetry_number, is_linear, E, S, Cv)
      !! Compute rotational contributions to thermodynamic properties.
      !!
      !! Uses rigid rotor approximation (classical limit, high T).
      !! For nonlinear: E = 3/2 RT, Cv = 3/2 R
      !! For linear: E = RT, Cv = R
      real(dp), intent(in) :: moments(3)      !! Principal moments in amu*Angstrom^2
      real(dp), intent(in) :: temperature     !! Temperature in K
      integer, intent(in) :: symmetry_number  !! Rotational symmetry number
      logical, intent(in) :: is_linear        !! True if linear molecule
      real(dp), intent(out) :: E              !! Thermal energy in Hartree
      real(dp), intent(out) :: S              !! Entropy in cal/(mol*K)
      real(dp), intent(out) :: Cv             !! Heat capacity in cal/(mol*K)

      real(dp) :: theta_rot(3)
      real(dp) :: T, sigma
      real(dp) :: qr
      integer :: i

      T = temperature
      sigma = real(symmetry_number, dp)

      ! Calculate rotational temperatures: theta_rot = h^2 / (8*pi^2*I*k_B)
      ! For I in amu*Angstrom^2:
      !   I_SI = I * 1.66054e-47 kg*m^2
      !   theta_rot = h^2 / (8*pi^2 * I_SI * k_B)
      !             = (6.62607e-34)^2 / (8 * pi^2 * 1.66054e-47 * 1.38065e-23 * I)
      !             = ROTTEMP_AMUA2_TO_K / I  [K]
      do i = 1, 3
         if (moments(i) > LINEAR_THRESHOLD) then
            theta_rot(i) = ROTTEMP_AMUA2_TO_K/moments(i)
         else
            theta_rot(i) = 0.0_dp
         end if
      end do

      if (is_linear) then
         ! Linear molecule: 2 rotational degrees of freedom
         E = R_HARTREE*T
         Cv = R_CALMOLK

         ! S = R * [1 + ln(T / (sigma * theta_rot))]
         ! Use the average of the two non-zero theta values
         if (theta_rot(3) > 0.0_dp) then
            qr = T/(sigma*theta_rot(3))
            S = R_CALMOLK*(1.0_dp + log(qr))
         else
            S = 0.0_dp
         end if
      else
         ! Nonlinear molecule: 3 rotational degrees of freedom
         E = 1.5_dp*R_HARTREE*T
         Cv = 1.5_dp*R_CALMOLK

         ! S = R * [3/2 + ln(sqrt(pi) * T^(3/2) / (sigma * sqrt(theta_A*theta_B*theta_C)))]
         if (theta_rot(1) > 0.0_dp .and. theta_rot(2) > 0.0_dp .and. theta_rot(3) > 0.0_dp) then
            qr = sqrt(PI)*(T**1.5_dp)/(sigma*sqrt(theta_rot(1)*theta_rot(2)*theta_rot(3)))
            S = R_CALMOLK*(1.5_dp + log(qr))
         else
            S = 0.0_dp
         end if
      end if

   end subroutine compute_rotational_thermo

   subroutine compute_vibrational_thermo(frequencies, n_freqs, temperature, E, S, Cv)
      !! Compute vibrational contributions to thermodynamic properties.
      !!
      !! Uses quantum harmonic oscillator partition function.
      !! E_vib = R * sum(theta_v * [1/(exp(u)-1)]) where u = theta_v/T
      !! S_vib = R * sum([u/(exp(u)-1) - ln(1-exp(-u))])
      !! Cv_vib = R * sum(u^2 * exp(u) / (exp(u)-1)^2)
      !!
      !! Note: ZPE is NOT included here (computed separately).
      real(dp), intent(in) :: frequencies(:)  !! Frequencies in cm^-1
      integer, intent(in) :: n_freqs          !! Number of frequencies
      real(dp), intent(in) :: temperature     !! Temperature in K
      real(dp), intent(out) :: E              !! Thermal energy in Hartree (excluding ZPE)
      real(dp), intent(out) :: S              !! Entropy in cal/(mol*K)
      real(dp), intent(out) :: Cv             !! Heat capacity in cal/(mol*K)

      integer :: i
      real(dp) :: T, freq, theta_v, u
      real(dp) :: exp_u, exp_neg_u
      real(dp) :: E_sum, S_sum, Cv_sum

      T = temperature
      E_sum = 0.0_dp
      S_sum = 0.0_dp
      Cv_sum = 0.0_dp

      do i = 1, n_freqs
         freq = frequencies(i)

         ! Skip imaginary and near-zero frequencies
         if (freq <= IMAG_FREQ_THRESHOLD) cycle
         if (freq < 10.0_dp) cycle  ! Skip very low frequencies (likely trans/rot residuals)

         ! Vibrational temperature: theta_v = h*c*nu / k = 1.4388 * nu (cm^-1)
         theta_v = CM1_TO_KELVIN*freq

         ! Reduced temperature ratio
         u = theta_v/T

         ! Avoid numerical issues for very large u (very low T or high freq)
         if (u > VIB_CLASSICAL_LIMIT) then
            ! Classical limit: modes are frozen out
            cycle
         end if

         exp_u = exp(u)
         exp_neg_u = exp(-u)

         ! Energy contribution (excluding ZPE): theta_v / (exp(u) - 1)
         E_sum = E_sum + theta_v/(exp_u - 1.0_dp)

         ! Entropy contribution: u/(exp(u)-1) - ln(1-exp(-u))
         S_sum = S_sum + u/(exp_u - 1.0_dp) - log(1.0_dp - exp_neg_u)

         ! Heat capacity contribution: u^2 * exp(u) / (exp(u)-1)^2
         Cv_sum = Cv_sum + (u*u*exp_u)/((exp_u - 1.0_dp)**2)
      end do

      ! Convert to proper units
      E = KB_HARTREE*E_sum           ! Hartree
      S = R_CALMOLK*S_sum            ! cal/(mol*K)
      Cv = R_CALMOLK*Cv_sum          ! cal/(mol*K)

   end subroutine compute_vibrational_thermo

   subroutine compute_electronic_entropy(spin_multiplicity, S_elec)
      !! Compute electronic entropy contribution.
      !!
      !! S_elec = R * ln(2S+1) where 2S+1 is the spin multiplicity.
      !! For singlet (mult=1): S_elec = 0
      integer, intent(in) :: spin_multiplicity  !! Electronic spin multiplicity (2S+1)
      real(dp), intent(out) :: S_elec           !! Electronic entropy in cal/(mol*K)

      S_elec = R_CALMOLK*log(real(spin_multiplicity, dp))

   end subroutine compute_electronic_entropy

   subroutine compute_partition_functions(total_mass, moments, frequencies, n_freqs, &
                                          temperature, pressure, sigma, is_linear, &
                                          q_trans, q_rot, q_vib)
      !! Compute partition functions for translation, rotation, and vibration.
      real(dp), intent(in) :: total_mass        !! Total mass in amu
      real(dp), intent(in) :: moments(3)        !! Principal moments in amu*Angstrom^2
      real(dp), intent(in) :: frequencies(:)    !! Frequencies in cm^-1
      integer, intent(in) :: n_freqs            !! Number of frequencies
      real(dp), intent(in) :: temperature       !! Temperature in K
      real(dp), intent(in) :: pressure          !! Pressure in atm
      integer, intent(in) :: sigma              !! Symmetry number
      logical, intent(in) :: is_linear          !! True if linear molecule
      real(dp), intent(out) :: q_trans          !! Translational partition function
      real(dp), intent(out) :: q_rot            !! Rotational partition function
      real(dp), intent(out) :: q_vib            !! Vibrational partition function

      real(dp) :: mass_kg, T, P_pa
      real(dp) :: lambda, V_molar
      real(dp) :: theta_rot(3), u
      integer :: i

      T = temperature
      mass_kg = total_mass*AMU_TO_KG
      P_pa = pressure*ATM_TO_PA

      ! Translational partition function: q_trans = (2*pi*m*k*T/h^2)^(3/2) * V
      ! where V = kT/P for ideal gas (per molecule)
      lambda = H_SI/sqrt(2.0_dp*PI*mass_kg*KB_SI*T)  ! thermal de Broglie wavelength
      V_molar = KB_SI*T/P_pa  ! volume per molecule
      q_trans = V_molar/(lambda**3)

      ! Rotational partition function
      ! theta_rot = h^2 / (8*pi^2*I*k_B) = ROTTEMP_AMUA2_TO_K / I (for I in amu*Angstrom^2)
      do i = 1, 3
         if (moments(i) > 1.0e-6_dp) then
            theta_rot(i) = ROTTEMP_AMUA2_TO_K/moments(i)
         else
            theta_rot(i) = 0.0_dp
         end if
      end do

      if (is_linear) then
         ! Linear: q_rot = T / (sigma * theta_rot)
         if (theta_rot(3) > 0.0_dp) then
            q_rot = T/(real(sigma, dp)*theta_rot(3))
         else
            q_rot = 1.0_dp
         end if
      else
         ! Nonlinear: q_rot = sqrt(pi) * T^(3/2) / (sigma * sqrt(theta_A * theta_B * theta_C))
         if (theta_rot(1) > 0.0_dp .and. theta_rot(2) > 0.0_dp .and. theta_rot(3) > 0.0_dp) then
            q_rot = sqrt(PI)*(T**1.5_dp)/ &
                    (real(sigma, dp)*sqrt(theta_rot(1)*theta_rot(2)*theta_rot(3)))
         else
            q_rot = 1.0_dp
         end if
      end if

      ! Vibrational partition function: q_vib = Product_i [1 / (1 - exp(-u_i))]
      ! where u_i = h*nu/(k*T) = theta_vib/T = 1.4388 * nu(cm^-1) / T
      q_vib = 1.0_dp
      do i = 1, n_freqs
         if (frequencies(i) > 10.0_dp) then  ! Skip near-zero frequencies
            u = CM1_TO_KELVIN*frequencies(i)/T
            if (u < 100.0_dp) then  ! Avoid overflow
               q_vib = q_vib/(1.0_dp - exp(-u))
            end if
         end if
      end do

   end subroutine compute_partition_functions

   subroutine compute_thermochemistry(coords, atomic_numbers, frequencies, n_atoms, n_freqs, &
                                      result, temperature, pressure, symmetry_number, spin_multiplicity)
      !! Main driver for thermochemistry calculations.
      !!
      !! Computes all thermodynamic quantities from molecular geometry and vibrational frequencies.
      real(dp), intent(in) :: coords(:, :)           !! Coordinates (3, n_atoms) in Bohr
      integer, intent(in) :: atomic_numbers(:)       !! Atomic numbers
      real(dp), intent(in) :: frequencies(:)         !! Frequencies in cm^-1
      integer, intent(in) :: n_atoms                 !! Number of atoms
      integer, intent(in) :: n_freqs                 !! Number of frequencies
      type(thermochemistry_result_t), intent(out) :: result  !! Output results
      real(dp), intent(in), optional :: temperature  !! Temperature in K (default 298.15)
      real(dp), intent(in), optional :: pressure     !! Pressure in atm (default 1.0)
      integer, intent(in), optional :: symmetry_number    !! Symmetry number (default 1)
      integer, intent(in), optional :: spin_multiplicity  !! Spin multiplicity (default 1)

      real(dp) :: center_of_mass(3)
      real(dp) :: principal_axes(3, 3)
      real(dp) :: T, P
      integer :: sigma, mult
      real(dp) :: S_total

      ! Set parameters
      T = DEFAULT_TEMPERATURE
      P = DEFAULT_PRESSURE
      sigma = DEFAULT_SYMMETRY_NUMBER
      mult = DEFAULT_SPIN_MULTIPLICITY

      if (present(temperature)) T = temperature
      if (present(pressure)) P = pressure
      if (present(symmetry_number)) sigma = symmetry_number
      if (present(spin_multiplicity)) mult = spin_multiplicity

      result%temperature = T
      result%pressure = P
      result%symmetry_number = sigma
      result%spin_multiplicity = mult

      ! Compute moments of inertia
      call compute_moments_of_inertia(coords, atomic_numbers, n_atoms, &
                                      center_of_mass, result%moments, principal_axes, &
                                      result%is_linear, result%total_mass)

      ! Compute rotational constants
      call compute_rotational_constants(result%moments, result%is_linear, result%rot_const)

      ! Compute ZPE
      call compute_zpe(frequencies, n_freqs, result%n_real_freqs, &
                       result%zpe_hartree, result%zpe_kcalmol)
      result%n_imag_freqs = n_freqs - result%n_real_freqs

      ! Compute translational contributions
      call compute_translational_thermo(result%total_mass, T, P, &
                                        result%E_trans, result%S_trans, result%Cv_trans)

      ! Compute rotational contributions
      call compute_rotational_thermo(result%moments, T, sigma, result%is_linear, &
                                     result%E_rot, result%S_rot, result%Cv_rot)

      ! Compute vibrational contributions (thermal, excluding ZPE)
      call compute_vibrational_thermo(frequencies, n_freqs, T, &
                                      result%E_vib, result%S_vib, result%Cv_vib)

      ! Compute electronic contributions
      result%E_elec = 0.0_dp  ! Ground state only
      call compute_electronic_entropy(mult, result%S_elec)

      ! Compute partition functions
      call compute_partition_functions(result%total_mass, result%moments, frequencies, n_freqs, &
                                       T, P, sigma, result%is_linear, &
                                       result%q_trans, result%q_rot, result%q_vib)

      ! Compute summary quantities
      ! Thermal correction to energy = ZPE + E_trans + E_rot + E_vib + E_elec
      result%thermal_correction_energy = result%zpe_hartree + &
                                         result%E_trans + result%E_rot + result%E_vib + result%E_elec

      ! Thermal correction to enthalpy = E_tot + RT
      result%thermal_correction_enthalpy = result%thermal_correction_energy + R_HARTREE*T

      ! Total entropy in Hartree/K
      S_total = (result%S_trans + result%S_rot + result%S_vib + result%S_elec)/HARTREE_TO_CALMOL

      ! Thermal correction to Gibbs = H - TS
      result%thermal_correction_gibbs = result%thermal_correction_enthalpy - T*S_total

   end subroutine compute_thermochemistry

   subroutine print_thermochemistry(result, electronic_energy, unit)
      !! Print thermochemistry results.
      type(thermochemistry_result_t), intent(in) :: result
      real(dp), intent(in) :: electronic_energy  !! Electronic energy in Hartree
      integer, intent(in), optional :: unit      !! Output unit (default: stdout)

      real(dp) :: H_total_cal, Cv_total, S_total, S_total_J
      real(dp) :: H_vib_cal, H_rot_cal, H_trans_cal
      real(dp) :: H_T, TS, G_T
      real(dp) :: total_energy, total_enthalpy, total_free_energy
      real(dp) :: H0_HT_PV, pV_cal
      character(len=512) :: line

      ! pV term = RT for ideal gas (in cal/mol)
      pV_cal = R_CALMOLK*result%temperature

      ! Compute enthalpy contributions (what xtb shows in "enthalpy" column)
      ! VIB: thermal vibrational energy only (ZPE is reported separately)
      H_vib_cal = result%E_vib*HARTREE_TO_CALMOL
      ! ROT: just internal energy (3/2)RT for nonlinear, RT for linear - no pV term
      H_rot_cal = result%E_rot*HARTREE_TO_CALMOL
      ! TR: includes pV term, so (5/2)RT for translation enthalpy
      H_trans_cal = result%E_trans*HARTREE_TO_CALMOL + pV_cal
      ! TOT: sum of thermal contributions (VIB + ROT + TR), NOT including ZPE
      H_total_cal = H_vib_cal + H_rot_cal + H_trans_cal

      ! Use Cp for translation (Cv + R) as is standard for thermochemistry output
      Cv_total = (result%Cv_trans + R_CALMOLK) + result%Cv_rot + result%Cv_vib
      S_total = result%S_trans + result%S_rot + result%S_vib + result%S_elec
      S_total_J = S_total*CAL_TO_J  ! cal to J

      ! Thermodynamic quantities in Hartree
      ! H(0)-H(T)+PV is the thermal correction WITHOUT ZPE (just E_trans + E_rot + E_vib + RT)
      H0_HT_PV = result%E_trans + result%E_rot + result%E_vib + R_HARTREE*result%temperature
      ! H(T) is full thermal correction including ZPE
      H_T = result%thermal_correction_enthalpy
      TS = result%thermal_correction_enthalpy - result%thermal_correction_gibbs
      G_T = result%thermal_correction_gibbs

      total_energy = electronic_energy
      total_enthalpy = electronic_energy + result%thermal_correction_enthalpy
      total_free_energy = electronic_energy + result%thermal_correction_gibbs

      ! Print header
      call logger%info(" ")
      call logger%info("Thermochemistry (RRHO)")
      call logger%info("======================")
      call logger%info(" ")

      ! Setup section - simple list
      write (line, '(A,F10.4,A)') "  Temperature:       ", result%temperature, " K"
      call logger%info(trim(line))
      write (line, '(A,F10.4,A)') "  Pressure:          ", result%pressure, " atm"
      call logger%info(trim(line))
      write (line, '(A,F10.4,A)') "  Molecular mass:    ", result%total_mass, " amu"
      call logger%info(trim(line))
      write (line, '(A,I6)') "  Vibrational modes: ", result%n_real_freqs
      call logger%info(trim(line))
      if (result%n_imag_freqs > 0) then
         write (line, '(A,I6,A)') "  Imaginary freqs:   ", result%n_imag_freqs, " (skipped)"
         call logger%info(trim(line))
      end if
      if (result%is_linear) then
         call logger%info("  Linear molecule:   yes")
      else
         call logger%info("  Linear molecule:   no")
      end if
      write (line, '(A,I6)') "  Symmetry number:   ", result%symmetry_number
      call logger%info(trim(line))
      call logger%info(" ")

      ! Contribution table
      call logger%info("  temp (K)       q        H(cal/mol)  Cp(cal/K/mol)  S(cal/K/mol)  S(J/K/mol)")
      call logger%info("  -------------------------------------------------------------------------")

      write (line, '(F8.2,A,ES10.3,F12.3,F14.3,F14.3,F12.3)') &
         result%temperature, "  VIB", result%q_vib, H_vib_cal, result%Cv_vib, &
         result%S_vib, result%S_vib*CAL_TO_J
      call logger%info(trim(line))

      write (line, '(A,ES10.3,F12.3,F14.3,F14.3,F12.3)') &
         "          ROT", result%q_rot, H_rot_cal, result%Cv_rot, &
         result%S_rot, result%S_rot*CAL_TO_J
      call logger%info(trim(line))

      write (line, '(A,ES10.3,F12.3,F14.3,F14.3,F12.3)') &
         "          INT", result%q_rot*result%q_vib, H_vib_cal + H_rot_cal, &
         result%Cv_vib + result%Cv_rot, result%S_vib + result%S_rot, &
         (result%S_vib + result%S_rot)*CAL_TO_J
      call logger%info(trim(line))

      ! For TR, report Cp = Cv + R (constant pressure heat capacity for ideal gas)
      write (line, '(A,ES10.3,F12.3,F14.3,F14.3,F12.3)') &
         "          TR ", result%q_trans, H_trans_cal, result%Cv_trans + R_CALMOLK, &
         result%S_trans, result%S_trans*CAL_TO_J
      call logger%info(trim(line))

      call logger%info("  -------------------------------------------------------------------------")
      write (line, '(A,F12.3,F14.3,F14.3,F12.3)') &
         "          TOT           ", H_total_cal, Cv_total, S_total, S_total_J
      call logger%info(trim(line))

      call logger%info(" ")

      ! Thermal corrections table
      call logger%info(" ")
      call logger%info("Thermal Corrections (Hartree)")
      call logger%info("-----------------------------")
      write (line, '(A,F18.12)') "  Zero-point energy:              ", result%zpe_hartree
      call logger%info(trim(line))
      write (line, '(A,F18.12)') "  Thermal correction to Energy:   ", result%thermal_correction_energy
      call logger%info(trim(line))
      write (line, '(A,F18.12)') "  Thermal correction to Enthalpy: ", result%thermal_correction_enthalpy
      call logger%info(trim(line))
      write (line, '(A,F18.12)') "  Thermal correction to Gibbs:    ", result%thermal_correction_gibbs
      call logger%info(trim(line))
      call logger%info(" ")

      ! Final totals
      call logger%info("Total Energies (Hartree)")
      call logger%info("------------------------")
      write (line, '(A,F20.12)') "  Electronic energy:            ", total_energy
      call logger%info(trim(line))
      write (line, '(A,F20.12)') "  Electronic + ZPE:             ", total_energy + result%zpe_hartree
      call logger%info(trim(line))
      write (line, '(A,F20.12)') "  Electronic + thermal (E):     ", total_energy + result%thermal_correction_energy
      call logger%info(trim(line))
      write (line, '(A,F20.12)') "  Electronic + thermal (H):     ", total_enthalpy
      call logger%info(trim(line))
      write (line, '(A,F20.12)') "  Electronic + thermal (G):     ", total_free_energy
      call logger%info(trim(line))
      call logger%info(" ")

   end subroutine print_thermochemistry

end module mqc_thermochemistry
