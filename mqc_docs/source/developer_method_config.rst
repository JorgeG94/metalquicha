.. _developer_method_config:

======================================
Method Configuration Architecture
======================================

This document explains the design philosophy behind metalquicha's method
configuration system and provides guidance for developers who want to add
new quantum chemistry methods.

Design Philosophy
=================

The method configuration system follows several key principles:

1. **Initialize once, use everywhere**: All configuration is read from the
   input file at startup and stored in a single ``method_config_t`` structure.
   This eliminates scattered configuration and makes the data flow explicit.

2. **Composition over inheritance**: Rather than complex class hierarchies,
   we use composition with nested configuration types. Shared settings live
   at the top level or in dedicated shared types.

3. **Shared configurations for related methods**: Methods that share common
   parameters (e.g., HF and DFT both need SCF settings) use shared config
   types to avoid duplication and ensure consistency.

4. **Factory pattern for instantiation**: A single factory creates all method
   instances, making it easy to add new methods without touching calling code.

Configuration Type Hierarchy
============================

The main configuration type uses composition to organize settings::

    type :: method_config_t
       !----- Common settings (all methods) -----
       integer(int32) :: method_type    ! METHOD_TYPE_MP2, etc.
       character(len=32) :: basis_set   ! "cc-pvdz", "aug-cc-pvtz", etc.
       logical :: use_spherical         ! Spherical vs Cartesian gaussians
       logical :: verbose               ! Enable verbose output

       !----- Shared configurations -----
       type(scf_config_t) :: scf        ! Used by HF and DFT
       type(correlation_config_t) :: corr  ! Used by MP2, CC, etc.

       !----- Method-specific configurations -----
       type(xtb_config_t) :: xtb        ! Semi-empirical xTB
       type(dft_config_t) :: dft        ! DFT-specific (functional, grid)
       type(mcscf_config_t) :: mcscf    ! Multi-reference
       type(cc_config_t) :: cc          ! Coupled-cluster specific
       type(f12_config_t) :: f12        ! F12 explicitly correlated
    end type

Shared Configuration Types
--------------------------

**scf_config_t** - Shared by HF and DFT::

    type :: scf_config_t
       integer :: max_iter = 100
       real(dp) :: energy_convergence = 1.0e-8_dp
       real(dp) :: density_convergence = 1.0e-6_dp
       logical :: use_diis = .true.
       integer :: diis_size = 8
    end type

**correlation_config_t** - Shared by all post-HF methods (MP2, CC, etc.)::

    type :: correlation_config_t
       real(dp) :: energy_convergence = 1.0e-8_dp
       integer :: n_frozen_core = -1     ! -1 = auto from elements
       logical :: freeze_core = .true.
       logical :: use_df = .true.        ! Density fitting (RI)
       character(len=32) :: aux_basis    ! RI auxiliary basis
       logical :: use_local = .false.    ! Local correlation (DLPNO, etc.)
       character(len=16) :: local_type   ! "dlpno", "pno", "lmp2"
       real(dp) :: pno_threshold         ! PNO truncation threshold
       logical :: use_scs = .false.      ! Spin-component scaling
       real(dp) :: scs_ss, scs_os        ! SCS scaling factors
    end type

Method-Specific Configuration Types
-----------------------------------

**cc_config_t** - Coupled-cluster specific settings::

    type :: cc_config_t
       integer :: max_iter = 100
       real(dp) :: amplitude_convergence = 1.0e-7_dp
       logical :: include_triples = .false.   ! (T) correction
       logical :: perturbative_triples = .true.
       logical :: use_diis = .true.
       integer :: diis_size = 8
       integer :: n_roots = 0            ! EOM-CC roots (0 = ground only)
       character(len=8) :: eom_type      ! "ee", "ip", "ea"
    end type

**f12_config_t** - Explicitly correlated F12 settings::

    type :: f12_config_t
       real(dp) :: geminal_exponent = 1.0_dp
       character(len=8) :: ansatz = '3c'      ! "3c", "3c(fix)", "2b"
       character(len=32) :: cabs_basis        ! CABS auxiliary basis
       character(len=32) :: optri_basis       ! Optional RI basis
       logical :: use_exponent_fit = .false.
       logical :: scale_triples = .true.
    end type

Method Type Constants
=====================

Method types are defined as integer constants in ``mqc_method_types.f90``
with a logical numbering scheme that groups related methods::

    ! Semi-empirical (1-9)
    METHOD_TYPE_GFN1 = 1
    METHOD_TYPE_GFN2 = 2

    ! SCF methods (10-19)
    METHOD_TYPE_HF = 10
    METHOD_TYPE_DFT = 11

    ! Multi-reference (20-29)
    METHOD_TYPE_MCSCF = 20

    ! Perturbation theory (30-39)
    METHOD_TYPE_MP2 = 30
    METHOD_TYPE_MP2_F12 = 31

    ! Coupled cluster (40-59)
    METHOD_TYPE_CCSD = 40
    METHOD_TYPE_CCSD_T = 41
    METHOD_TYPE_CCSD_F12 = 42
    METHOD_TYPE_CCSD_T_F12 = 43

This numbering scheme leaves room for future methods (e.g., MP3 at 32,
CC2 at 44, CCSDT at 45).

The Factory Pattern
===================

The ``method_factory_t`` in ``mqc_method_factory.F90`` creates method
instances from configuration::

    function factory_create(this, config) result(method)
       class(method_factory_t), intent(in) :: this
       type(method_config_t), intent(in) :: config
       class(qc_method_t), allocatable :: method

       select case (config%method_type)
       case (METHOD_TYPE_HF)
          allocate(hf_method_t :: method)
          call configure_hf(method, config)

       case (METHOD_TYPE_MP2)
          allocate(mp2_method_t :: method)
          call configure_mp2(method, config)

       case (METHOD_TYPE_CCSD, METHOD_TYPE_CCSD_T)
          allocate(cc_method_t :: method)
          call configure_cc(method, config)
       ! ... etc
       end select
    end function

Each ``configure_*`` subroutine reads from the appropriate config sections::

    subroutine configure_cc(method, config)
       ! Common settings
       m%basis_set = config%basis_set

       ! SCF settings (for reference calculation)
       m%scf_max_iter = config%scf%max_iter

       ! Correlation settings (shared)
       m%freeze_core = config%corr%freeze_core
       m%use_df = config%corr%use_df
       m%aux_basis = config%corr%aux_basis

       ! CC-specific settings
       m%max_iter = config%cc%max_iter
       m%include_triples = config%cc%include_triples
    end subroutine

Adding a New Method
===================

To add a new quantum chemistry method, follow these steps:

Step 1: Add Method Type Constant
--------------------------------

In ``src/mqc_method_types.f90``:

1. Add the constant with an appropriate number::

    integer(int32), parameter :: METHOD_TYPE_MP3 = 32

2. Export it in the public statement::

    public :: METHOD_TYPE_MP3

3. Add string conversions in ``method_type_from_string`` and
   ``method_type_to_string``::

    case ('mp3', 'ri-mp3')
       method_type = METHOD_TYPE_MP3

    case (METHOD_TYPE_MP3)
       method_str = "mp3"

Step 2: Add Configuration (if needed)
-------------------------------------

If your method needs settings not covered by existing config types:

1. **Prefer using existing shared configs** when possible. For example,
   MP3 would just use ``correlation_config_t`` with no additions.

2. **Add to an existing config type** if the setting is broadly applicable.
   For example, adding a ``use_t1_diagnostic`` flag to ``cc_config_t``.

3. **Create a new config type** only if the method has truly unique settings.
   Add it to ``mqc_method_config.f90``::

    type :: mp3_config_t
       logical :: use_laplace = .false.  ! Laplace transform MP3
       integer :: laplace_points = 5
    end type

   Then add it to ``method_config_t`` and ``config_reset``.

Step 3: Create the Method Type
------------------------------

Create ``src/methods/mqc_method_mp3.f90``::

    module mqc_method_mp3
       use mqc_method_base, only: qc_method_t
       implicit none

       type :: mp3_options_t
          ! Copy relevant settings from config
          character(len=32) :: basis_set
          logical :: freeze_core
          ! ... etc
       end type

       type, extends(qc_method_t) :: mp3_method_t
          type(mp3_options_t) :: options
       contains
          procedure :: calc_energy => mp3_calc_energy
          procedure :: calc_gradient => mp3_calc_gradient
          procedure :: calc_hessian => mp3_calc_hessian
       end type

    contains
       ! Implement the methods...
    end module

Step 4: Add to the Factory
--------------------------

In ``src/methods/mqc_method_factory.F90``:

1. Add the use statement::

    use mqc_method_mp3, only: mp3_method_t

2. Add the case in ``factory_create``::

    case (METHOD_TYPE_MP3)
       allocate(mp3_method_t :: method)
       call configure_mp3(method, config)

3. Implement ``configure_mp3``::

    subroutine configure_mp3(method, config)
       class(qc_method_t), intent(inout) :: method
       type(method_config_t), intent(in) :: config

       select type (m => method)
       type is (mp3_method_t)
          m%options%basis_set = config%basis_set
          m%options%freeze_core = config%corr%freeze_core
          m%options%use_df = config%corr%use_df
          ! ... etc
       end select
    end subroutine

Step 5: Add to CMakeLists.txt
-----------------------------

In ``src/methods/CMakeLists.txt``::

    target_sources(${main_lib} PRIVATE mqc_method_mp3.f90)

Step 6: Update Documentation
----------------------------

Update this document and any user-facing documentation to include the
new method.

Example: Method Configuration in Practice
=========================================

Here's how different calculation types use the configuration:

**Simple MP2**::

    config%method_type = METHOD_TYPE_MP2
    config%basis_set = 'cc-pvtz'
    config%corr%freeze_core = .true.
    config%corr%use_df = .true.
    config%corr%aux_basis = 'cc-pvtz-ri'

**DLPNO-CCSD(T)**::

    config%method_type = METHOD_TYPE_CCSD_T
    config%basis_set = 'cc-pvtz'
    config%corr%freeze_core = .true.
    config%corr%use_local = .true.
    config%corr%local_type = 'dlpno'
    config%corr%pno_threshold = 1.0e-7_dp
    config%cc%include_triples = .true.

**CCSD(T)-F12**::

    config%method_type = METHOD_TYPE_CCSD_T_F12
    config%basis_set = 'cc-pvtz-f12'
    config%corr%freeze_core = .true.
    config%cc%include_triples = .true.
    config%f12%ansatz = '3c'
    config%f12%geminal_exponent = 1.0_dp
    config%f12%cabs_basis = 'cc-pvtz-f12-cabs'

**EOM-CCSD for excited states**::

    config%method_type = METHOD_TYPE_CCSD
    config%basis_set = 'aug-cc-pvdz'
    config%cc%n_roots = 5
    config%cc%eom_type = 'ee'

Design Rationale
================

Why Composition Over Inheritance?
---------------------------------

Fortran's OOP support for inheritance is functional but can be awkward,
especially for configuration types. Composition provides:

- **Explicit data flow**: It's clear which settings come from where
- **Flexibility**: A method can use any combination of shared configs
- **Simplicity**: No virtual dispatch overhead for config access
- **Extensibility**: Adding new config types doesn't affect existing ones

Why Shared Config Types?
------------------------

Many quantum chemistry methods share common concepts:

- HF and DFT both solve SCF equations -> ``scf_config_t``
- MP2 and CCSD both have frozen core, RI, local approximations -> ``correlation_config_t``
- All F12 methods share geminal/CABS settings -> ``f12_config_t``

Sharing these prevents duplication and ensures that common settings are
handled consistently across methods.

Why a Factory?
--------------

The factory pattern provides:

- **Single point of method creation**: Easy to add new methods
- **Encapsulated instantiation**: Calling code doesn't need to know
  concrete method types
- **Consistent configuration**: All methods are configured the same way
- **Testability**: Easy to mock or substitute methods

Why Integer Constants for Method Types?
---------------------------------------

Integer constants (vs. strings) provide:

- **Fast comparison**: No string operations in hot paths
- **Compile-time checking**: Typos are caught by the compiler
- **Clear grouping**: The numbering scheme documents method families
- **Memory efficiency**: 4 bytes vs. variable-length strings

File Locations
==============

- ``src/mqc_method_types.f90`` - Method type constants and conversions
- ``src/methods/mqc_method_config.f90`` - All configuration types
- ``src/methods/mqc_method_factory.F90`` - Factory implementation
- ``src/methods/mqc_method_base.f90`` - Abstract base class
- ``src/methods/mqc_method_*.f90`` - Concrete method implementations
