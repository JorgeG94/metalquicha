!! Backend-independent description of a cuEST SCF calculation
module mqc_cuest_iface
   !! Holds `cuest_scf_settings_t`, the settings an HF or DFT method hands to
   !! the cuEST backend.
   !!
   !! It lives here, in `src/`, rather than with the backend for one reason:
   !! the backend cannot be part of the fpm build (fpm cannot link the cuEST
   !! library, and its module-dependency scanner is preprocessor-blind, so it
   !! would demand the backend modules even from behind an `#ifdef`). Keeping
   !! the *type* here and reaching the backend through `mqc_cuest_bridge` --
   !! which has a stub form fpm compiles and a real form CMake compiles -- lets
   !! the method files carry no preprocessor conditionals at all.
   use pic_types, only: dp
   implicit none
   private

   public :: cuest_scf_settings_t

   type :: cuest_scf_settings_t
      !! Method-independent description of one cuEST SCF calculation
      character(len=32) :: basis_set = "sto-3g"
         !! Orbital basis set name
      character(len=32) :: aux_basis_set = "def2-universal-jkfit"
      logical :: density_fitting = .false.
         !! Read by the CPU backend only; cuEST always fits.
         !! Auxiliary (JKFIT) basis. Required: cuEST fits J and K always.
      character(len=32) :: functional = ""
         !! Exchange-correlation functional; empty means Hartree-Fock
      logical :: spherical = .true.
         !! Pure (spherical) vs Cartesian angular functions
      logical :: verbose = .false.
         !! Print the SCF iteration table
      character(len=16) :: guess = "gwh"
         !! Initial guess: 'core', 'gwh' or 'sac'
      integer :: device_rank = 0
         !! Node-local MPI rank; decides which GPU this rank binds to
      logical :: unrestricted = .false.
         !! Force UHF/UKS even for a closed shell

      logical :: allow_crap_scf = .false.
         !! Keep a non-converged SCF instead of failing the fragment. Same
         !! meaning and same default as the xTB path -- one physical condition
         !! must not have two behaviours depending on which backend ran it.
      integer :: max_iter = 100
      real(dp) :: energy_tol = 1.0e-8_dp
      real(dp) :: density_tol = 1.0e-6_dp
      logical :: use_diis = .true.
      integer :: diis_size = 8

      integer :: radial_points = 75    !! XC grid radial points per atom
      integer :: angular_points = 302  !! XC grid Lebedev order
   end type cuest_scf_settings_t

end module mqc_cuest_iface
