! ============================================================================
!  cuest.f90 -- Fortran 2008 iso_c_binding interface to NVIDIA cuEST
!
!  GENERATED FILE -- do not edit by hand.
!  Regenerate with:  python3 generate_cuest_fortran.py
!  Source: the cuEST C headers (include/cuest.h, v0.2.0).
!
!  Conventions
!  -----------
!  * Every cuEST function returns cuestStatus_t; here each is an
!    INTEGER(c_int) FUNCTION whose result is the status code.
!  * All opaque handles (cuestHandle_t, plans, bases, parameter and results
!    objects, ...) are C  void*  and map to TYPE(c_ptr).
!      - handle passed IN  -> TYPE(c_ptr), VALUE
!      - handle returned   -> TYPE(c_ptr), INTENT(OUT)   (C void**)
!  * Bulk numeric buffers (double*, e.g. matrices/coordinates) are, in cuEST,
!    generally GPU DEVICE pointers.  They are exposed as TYPE(c_ptr), VALUE so
!    you may pass either a device address or C_LOC(host_array).
!  * The two workspace structs are interoperable derived types (below).
!  * Enumerators are PUBLIC INTEGER(c_int) PARAMETERs with their C names.
!
!  Build:  gfortran -c cuest.f90     (produces cuest.mod)
!  Link :  ... -L<pkg>/lib -lcuest -Wl,-rpath,<pkg>/lib
! ============================================================================
module cuest
   use, intrinsic :: iso_c_binding
   implicit none
   public

   ! ---- library version (compile-time, from the headers) ----------------
   integer(c_int), parameter :: CUEST_VER_MAJOR = 0
   integer(c_int), parameter :: CUEST_VER_MINOR = 2
   integer(c_int), parameter :: CUEST_VER_PATCH = 0

   ! ======================================================================
   !  Enumerations (cuestStatus_t, handle/attribute ids, modes, ...)
   ! ======================================================================
   ! ---- cuestStatus_t
   integer(c_int), parameter :: CUEST_STATUS_SUCCESS = 0
   integer(c_int), parameter :: CUEST_STATUS_EXCEPTION = 1
   integer(c_int), parameter :: CUEST_STATUS_NULL_POINTER = 2
   integer(c_int), parameter :: CUEST_STATUS_INVALID_ARGUMENT = 3
   integer(c_int), parameter :: CUEST_STATUS_INVALID_SIZE = 4
   integer(c_int), parameter :: CUEST_STATUS_INVALID_TYPE = 5
   integer(c_int), parameter :: CUEST_STATUS_INVALID_PARAMETER = 6
   integer(c_int), parameter :: CUEST_STATUS_INVALID_ATTRIBUTE = 7
   integer(c_int), parameter :: CUEST_STATUS_INVALID_HANDLE = 8
   integer(c_int), parameter :: CUEST_STATUS_UNKNOWN_ERROR = 9
   integer(c_int), parameter :: CUEST_STATUS_UNSUPPORTED_ARGUMENT = 10
   integer(c_int), parameter :: CUEST_STATUS_UNSUPPORTED_ARCHITECTURE = 11
   integer(c_int), parameter :: CUEST_STATUS_INVALID_PLAN = 12
   integer(c_int), parameter :: CUEST_STATUS_HOME_NOT_FOUND = 13

   ! ---- cuestType_t
   integer(c_int), parameter :: CUEST_HANDLE = 0
   integer(c_int), parameter :: CUEST_AOSHELL = 1
   integer(c_int), parameter :: CUEST_AOBASIS = 2
   integer(c_int), parameter :: CUEST_AOPAIRLIST = 3
   integer(c_int), parameter :: CUEST_OEINTPLAN = 4
   integer(c_int), parameter :: CUEST_DFINTPLAN = 5
   integer(c_int), parameter :: CUEST_ATOMGRID = 6
   integer(c_int), parameter :: CUEST_MOLECULARGRID = 7
   integer(c_int), parameter :: CUEST_XCINTPLAN = 8
   integer(c_int), parameter :: CUEST_PCMINTPLAN = 9
   integer(c_int), parameter :: CUEST_ECPSHELL = 10
   integer(c_int), parameter :: CUEST_ECPATOM = 11
   integer(c_int), parameter :: CUEST_ECPINTPLAN = 12

   ! ---- cuestHandleAttributes_t
   integer(c_int), parameter :: CUEST_HANDLE_NONE = 0

   ! ---- cuestAOShellAttributes_t
   integer(c_int), parameter :: CUEST_AOSHELL_IS_PURE = 0
   integer(c_int), parameter :: CUEST_AOSHELL_L = 1
   integer(c_int), parameter :: CUEST_AOSHELL_NUM_PRIMITIVE = 2
   integer(c_int), parameter :: CUEST_AOSHELL_NUM_AO = 3
   integer(c_int), parameter :: CUEST_AOSHELL_NUM_PURE = 4
   integer(c_int), parameter :: CUEST_AOSHELL_NUM_CART = 5

   ! ---- cuestAOBasisAttributes_t
   integer(c_int), parameter :: CUEST_AOBASIS_IS_PURE = 0
   integer(c_int), parameter :: CUEST_AOBASIS_IS_CART = 1
   integer(c_int), parameter :: CUEST_AOBASIS_IS_MIXED = 2
   integer(c_int), parameter :: CUEST_AOBASIS_NUM_ATOM = 3
   integer(c_int), parameter :: CUEST_AOBASIS_NUM_SHELL = 4
   integer(c_int), parameter :: CUEST_AOBASIS_NUM_AO = 5
   integer(c_int), parameter :: CUEST_AOBASIS_NUM_CART = 6
   integer(c_int), parameter :: CUEST_AOBASIS_NUM_PRIMITIVE = 7
   integer(c_int), parameter :: CUEST_AOBASIS_MAX_L = 8

   ! ---- cuestAOPairListAttributes_t
   integer(c_int), parameter :: CUEST_AOPAIRLIST_NONE = 0

   ! ---- cuestECPShellAttributes_t
   integer(c_int), parameter :: CUEST_ECPSHELL_L = 0
   integer(c_int), parameter :: CUEST_ECPSHELL_NUM_PRIMITIVE = 1

   ! ---- cuestECPAtomAttributes_t
   integer(c_int), parameter :: CUEST_ECPATOM_MAX_L = 0
   integer(c_int), parameter :: CUEST_ECPATOM_NUM_ELECTRON = 1

   ! ---- cuestOEIntPlanAttributes_t
   integer(c_int), parameter :: CUEST_OEINTPLAN_NONE = 0

   ! ---- cuestECPIntPlanAttributes_t
   integer(c_int), parameter :: CUEST_ECPINTPLAN_NUM_ACTIVE_ATOM = 0
   integer(c_int), parameter :: CUEST_ECPINTPLAN_NUM_ELECTRON = 1
   integer(c_int), parameter :: CUEST_ECPINTPLAN_MAX_L = 2
   integer(c_int), parameter :: CUEST_ECPINTPLAN_NUM_RADIAL_POINT = 3
   integer(c_int), parameter :: CUEST_ECPINTPLAN_NUM_SPHERICAL_POINT = 4
   integer(c_int), parameter :: CUEST_ECPINTPLAN_NUM_POINT = 5

   ! ---- cuestDFIntPlanAttributes_t
   integer(c_int), parameter :: CUEST_DFINTPLAN_NONE = 0

   ! ---- cuestAtomGridAttributes_t
   integer(c_int), parameter :: CUEST_ATOMGRID_NUM_POINT = 0
   integer(c_int), parameter :: CUEST_ATOMGRID_NUM_RADIAL_POINT = 1
   integer(c_int), parameter :: CUEST_ATOMGRID_MAX_ANGULAR_POINT = 2

   ! ---- cuestMolecularGridAttributes_t
   integer(c_int), parameter :: CUEST_MOLECULARGRID_NUM_ATOM = 0
   integer(c_int), parameter :: CUEST_MOLECULARGRID_NUM_POINT = 1
   integer(c_int), parameter :: CUEST_MOLECULARGRID_MAX_POINT = 2
   integer(c_int), parameter :: CUEST_MOLECULARGRID_NUM_RADIAL_POINT = 3
   integer(c_int), parameter :: CUEST_MOLECULARGRID_MAX_RADIAL_POINT = 4
   integer(c_int), parameter :: CUEST_MOLECULARGRID_MAX_ANGULAR_POINT = 5

   ! ---- cuestXCIntPlanAttributes_t
   integer(c_int), parameter :: CUEST_XCINTPLAN_EXCHANGE_SCALE = 0
   integer(c_int), parameter :: CUEST_XCINTPLAN_LRC_EXCHANGE_SCALE = 1
   integer(c_int), parameter :: CUEST_XCINTPLAN_LRC_OMEGA = 2
   integer(c_int), parameter :: CUEST_XCINTPLAN_VV10_SCALE = 3
   integer(c_int), parameter :: CUEST_XCINTPLAN_VV10_C = 4
   integer(c_int), parameter :: CUEST_XCINTPLAN_VV10_B = 5
   integer(c_int), parameter :: CUEST_XCINTPLAN_IS_HYBRID = 6
   integer(c_int), parameter :: CUEST_XCINTPLAN_IS_LRC_HYBRID = 7
   integer(c_int), parameter :: CUEST_XCINTPLAN_IS_VV10 = 8
   integer(c_int), parameter :: CUEST_XCINTPLAN_ENGINE_DESCRIPTION = 9
   integer(c_int), parameter :: CUEST_XCINTPLAN_ENGINE_CITATION = 10
   integer(c_int), parameter :: CUEST_XCINTPLAN_FUNCTIONAL_CITATION = 11
   integer(c_int), parameter :: CUEST_XCINTPLAN_FUNCTIONAL_DESCRIPTION = 12

   ! ---- cuestPCMIntPlanAttributes_t
   integer(c_int), parameter :: CUEST_PCMINTPLAN_NUM_POINT = 0
   integer(c_int), parameter :: CUEST_PCMINTPLAN_NUM_ACTIVE_POINT = 1

   ! ---- cuestParametersType_t
   integer(c_int), parameter :: CUEST_HANDLE_PARAMETERS = 0
   integer(c_int), parameter :: CUEST_AOSHELL_PARAMETERS = 1
   integer(c_int), parameter :: CUEST_AOBASIS_PARAMETERS = 2
   integer(c_int), parameter :: CUEST_AOPAIRLIST_PARAMETERS = 3
   integer(c_int), parameter :: CUEST_OEINTPLAN_PARAMETERS = 4
   integer(c_int), parameter :: CUEST_DFINTPLAN_PARAMETERS = 5
   integer(c_int), parameter :: CUEST_ATOMGRID_PARAMETERS = 6
   integer(c_int), parameter :: CUEST_MOLECULARGRID_PARAMETERS = 7
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS = 8
   integer(c_int), parameter :: CUEST_OVERLAPCOMPUTE_PARAMETERS = 9
   integer(c_int), parameter :: CUEST_OVERLAPDERIVATIVECOMPUTE_PARAMETERS = 10
   integer(c_int), parameter :: CUEST_KINETICCOMPUTE_PARAMETERS = 11
   integer(c_int), parameter :: CUEST_KINETICDERIVATIVECOMPUTE_PARAMETERS = 12
   integer(c_int), parameter :: CUEST_POTENTIALCOMPUTE_PARAMETERS = 13
   integer(c_int), parameter :: CUEST_POTENTIALDERIVATIVECOMPUTE_PARAMETERS = 14
   integer(c_int), parameter :: CUEST_DFCOULOMBCOMPUTE_PARAMETERS = 15
   integer(c_int), parameter :: CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS = 16
   integer(c_int), parameter :: CUEST_DFSYMMETRICDERIVATIVECOMPUTE_PARAMETERS = 17
   integer(c_int), parameter :: CUEST_XCPOTENTIALRKSCOMPUTE_PARAMETERS = 18
   integer(c_int), parameter :: CUEST_XCPOTENTIALUKSCOMPUTE_PARAMETERS = 19
   integer(c_int), parameter :: CUEST_XCDERIVATIVERKSCOMPUTE_PARAMETERS = 20
   integer(c_int), parameter :: CUEST_XCDERIVATIVEUKSCOMPUTE_PARAMETERS = 21
   integer(c_int), parameter :: CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS = 22
   integer(c_int), parameter :: CUEST_NONLOCALXCPOTENTIALUKSCOMPUTE_PARAMETERS = 23
   integer(c_int), parameter :: CUEST_NONLOCALXCDERIVATIVERKSCOMPUTE_PARAMETERS = 24
   integer(c_int), parameter :: CUEST_NONLOCALXCDERIVATIVEUKSCOMPUTE_PARAMETERS = 25
   integer(c_int), parameter :: CUEST_XCDENSITYCOMPUTE_PARAMETERS = 26
   integer(c_int), parameter :: CUEST_XCPOTENTIALCOMPUTE_PARAMETERS = 27
   integer(c_int), parameter :: CUEST_XCDERIVATIVECOMPUTE_PARAMETERS = 28
   integer(c_int), parameter :: CUEST_XCGRIDDERIVATIVECOMPUTE_PARAMETERS = 29
   integer(c_int), parameter :: CUEST_XCINTEGRATIONGRIDCOMPUTE_PARAMETERS = 30
   integer(c_int), parameter :: CUEST_XCINTEGRATIONWEIGHTCOMPUTE_PARAMETERS = 31
   integer(c_int), parameter :: CUEST_PCMINTPLAN_PARAMETERS = 32
   integer(c_int), parameter :: CUEST_PCMPOTENTIALCOMPUTE_PARAMETERS = 33
   integer(c_int), parameter :: CUEST_PCMDERIVATIVECOMPUTE_PARAMETERS = 34
   integer(c_int), parameter :: CUEST_ECPSHELL_PARAMETERS = 35
   integer(c_int), parameter :: CUEST_ECPATOM_PARAMETERS = 36
   integer(c_int), parameter :: CUEST_ECPINTPLAN_PARAMETERS = 37
   integer(c_int), parameter :: CUEST_ECPCOMPUTE_PARAMETERS = 38
   integer(c_int), parameter :: CUEST_ECPDERIVATIVECOMPUTE_PARAMETERS = 39
   integer(c_int), parameter :: CUEST_DFNONSYMMETRICEXCHANGECOMPUTE_PARAMETERS = 40
   integer(c_int), parameter :: CUEST_MULTIPOLECOMPUTE_PARAMETERS = 41
   integer(c_int), parameter :: CUEST_MULTIPOLEDERIVATIVECOMPUTE_PARAMETERS = 42
   integer(c_int), parameter :: CUEST_ANGULARMOMENTUMCOMPUTE_PARAMETERS = 43
   integer(c_int), parameter :: CUEST_ANGULARMOMENTUMDERIVATIVECOMPUTE_PARAMETERS = 44
   integer(c_int), parameter :: CUEST_NABLACOMPUTE_PARAMETERS = 45
   integer(c_int), parameter :: CUEST_NABLADERIVATIVECOMPUTE_PARAMETERS = 46
   integer(c_int), parameter :: CUEST_PCMINTEGRATIONGRIDCOMPUTE_PARAMETERS = 47
   integer(c_int), parameter :: CUEST_PCMINTEGRATIONWEIGHTCOMPUTE_PARAMETERS = 48
   integer(c_int), parameter :: CUEST_DFMOINTEGRALSCOMPUTE_PARAMETERS = 49
   integer(c_int), parameter :: CUEST_XCNONSYMMETRICDENSITYCOMPUTE_PARAMETERS = 50
   integer(c_int), parameter :: CUEST_XCNONSYMMETRICDERIVATIVECOMPUTE_PARAMETERS = 51
   integer(c_int), parameter :: CUEST_PCMRADIIDERIVATIVECOMPUTE_PARAMETERS = 52

   ! ---- cuestHandleParametersAttributes_t
   integer(c_int), parameter :: CUEST_HANDLE_PARAMETERS_CUDASTREAM = 0
   integer(c_int), parameter :: CUEST_HANDLE_PARAMETERS_CUBLAS = 1
   integer(c_int), parameter :: CUEST_HANDLE_PARAMETERS_CUSOLVER = 2
   integer(c_int), parameter :: CUEST_HANDLE_PARAMETERS_MAX_GAUSS_HERMITE = 3
   integer(c_int), parameter :: CUEST_HANDLE_PARAMETERS_MAX_L_SOLID_HARMONIC = 4
   integer(c_int), parameter :: CUEST_HANDLE_PARAMETERS_MAX_RYS = 5
   integer(c_int), parameter :: CUEST_HANDLE_PARAMETERS_RYS_SCHEME = 6
   integer(c_int), parameter :: CUEST_HANDLE_PARAMETERS_JIT_CACHE_DIR = 7
   integer(c_int), parameter :: CUEST_HANDLE_PARAMETERS_JIT_COMPILE_THREADS = 8

   ! ---- cuestHandleParametersRysScheme_t
   integer(c_int), parameter :: CUEST_HANDLE_PARAMETERS_RYS_SCHEME_06_21_2025 = 0

   ! ---- cuestMathMode_t
   integer(c_int), parameter :: CUEST_DEFAULT_MATH_MODE = 0
   integer(c_int), parameter :: CUEST_NATIVE_FP64_MATH_MODE = 1

   ! ---- cuestJITUsageMode_t
   integer(c_int), parameter :: CUEST_JIT_USAGE_MODE_ON = 0
   integer(c_int), parameter :: CUEST_JIT_USAGE_MODE_OFF = 1

   ! ---- cuestFfloatUsageMode_t
   integer(c_int), parameter :: CUEST_FFLOAT_USAGE_MODE_DEFAULT = 0
   integer(c_int), parameter :: CUEST_FFLOAT_USAGE_MODE_ON = 1
   integer(c_int), parameter :: CUEST_FFLOAT_USAGE_MODE_OFF = 2

   ! ---- cuestAOShellParametersAttributes_t
   integer(c_int), parameter :: CUEST_AOSHELL_PARAMETERS_NONE = 0

   ! ---- cuestAOBasisParametersAttributes_t
   integer(c_int), parameter :: CUEST_AOBASIS_PARAMETERS_NONE = 0

   ! ---- cuestAOPairListParametersAttributes_t
   integer(c_int), parameter :: CUEST_AOPAIRLIST_PARAMETERS_NONE = 0

   ! ---- cuestECPShellParametersAttributes_t
   integer(c_int), parameter :: CUEST_ECPSHELL_PARAMETERS_NONE = 0

   ! ---- cuestECPAtomParametersAttributes_t
   integer(c_int), parameter :: CUEST_ECPATOM_PARAMETERS_NONE = 0

   ! ---- cuestOEIntPlanParametersAttributes_t
   integer(c_int), parameter :: CUEST_OEINTPLAN_PARAMETERS_NONE = 0

   ! ---- cuestECPIntPlanParametersAttributes_t
   integer(c_int), parameter :: CUEST_ECPINTPLAN_PARAMETERS_THRESHOLD_COLLOCATION = 0
   integer(c_int), parameter :: CUEST_ECPINTPLAN_PARAMETERS_ECP_CUTOFF = 1
   integer(c_int), parameter :: CUEST_ECPINTPLAN_PARAMETERS_R0 = 2
   integer(c_int), parameter :: CUEST_ECPINTPLAN_PARAMETERS_MULTIPLIER = 3
   integer(c_int), parameter :: CUEST_ECPINTPLAN_PARAMETERS_RADIAL_NUM_POINT = 4
   integer(c_int), parameter :: CUEST_ECPINTPLAN_PARAMETERS_SPHERICAL_NUM_POINT = 5

   ! ---- cuestECPComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_ECPCOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestECPDerivativeComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_ECPDERIVATIVECOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestDFIntPlanParametersAttributes_t
   integer(c_int), parameter :: CUEST_DFINTPLAN_PARAMETERS_FITTING_CUTOFF = 0
   integer(c_int), parameter :: CUEST_DFINTPLAN_PARAMETERS_FITTING_RELATIVE_CONDITIONING = 1
   integer(c_int), parameter :: CUEST_DFINTPLAN_PARAMETERS_EXCHANGE_FRACTION = 2
   integer(c_int), parameter :: CUEST_DFINTPLAN_PARAMETERS_FITTING_ALGORITHM = 3
   integer(c_int), parameter :: CUEST_DFINTPLAN_PARAMETERS_LRC_EXCHANGE_FRACTION = 4
   integer(c_int), parameter :: CUEST_DFINTPLAN_PARAMETERS_LRC_EXCHANGE_OMEGA = 5
   integer(c_int), parameter :: CUEST_DFINTPLAN_PARAMETERS_THREE_INDEX_INTEGRAL_DIRECT = 6

   ! ---- cuestDFIntPlanParametersFittingAlgorithm_t
   integer(c_int), parameter :: CUEST_DFINTPLAN_PARAMETERS_FITTING_ALGORITHM_QR = 0
   integer(c_int), parameter :: CUEST_DFINTPLAN_PARAMETERS_FITTING_ALGORITHM_MATRIXPOWER = 1

   ! ---- cuestAtomGridParametersAttributes_t
   integer(c_int), parameter :: CUEST_ATOMGRID_PARAMETERS_NONE = 0

   ! ---- cuestMolecularGridParametersAttributes_t
   integer(c_int), parameter :: CUEST_MOLECULARGRID_PARAMETERS_NONE = 0

   ! ---- cuestXCIntPlanParametersAttributes_t
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_THRESHOLD_COLLOCATION = 0
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_DERIVATIVE_LEVEL = 1
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_NUM_THREADS = 2

   ! ---- cuestXCIntPlanParametersFunctional_t
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_HF = 0
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_B3LYP1 = 1
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_B3LYP5 = 2
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_B97 = 3
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_BLYP = 4
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_M06L = 5
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_PBE = 6
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_PBE0 = 7
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_R2SCAN = 8
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_SVWN5 = 9
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_B97MV = 10
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_LCWPBE = 11
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_WB97X = 12
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_WB97XV = 13
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_WB97MV = 14
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_LCWPBEH = 15
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_CAMB3LYP = 16
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_HSE06 = 17
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_M06 = 18
   integer(c_int), parameter :: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_M062X = 19

   ! ---- cuestOverlapComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_OVERLAPCOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestOverlapDerivativeComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_OVERLAPDERIVATIVECOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestKineticComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_KINETICCOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestKineticDerivativeComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_KINETICDERIVATIVECOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestPotentialComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_POTENTIALCOMPUTE_PARAMETERS_FFLOAT_USAGE_MODE = 0
   integer(c_int), parameter :: CUEST_POTENTIALCOMPUTE_PARAMETERS_JIT_USAGE_MODE = 1

   ! ---- cuestPotentialDerivativeComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_POTENTIALDERIVATIVECOMPUTE_PARAMETERS_FFLOAT_USAGE_MODE = 0
   integer(c_int), parameter :: CUEST_POTENTIALDERIVATIVECOMPUTE_PARAMETERS_JIT_USAGE_MODE = 1

   ! ---- cuestDFCoulombComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_DFCOULOMBCOMPUTE_PARAMETERS_FFLOAT_USAGE_MODE = 0
   integer(c_int), parameter :: CUEST_DFCOULOMBCOMPUTE_PARAMETERS_JIT_USAGE_MODE = 1

   ! ---- cuestDFSymmetricExchangeComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS_INT8_SLICE_COUNT = 0
   integer(c_int), parameter :: CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS_INT8_MODULUS_COUNT = 1

   ! ---- cuestDFNonsymmetricExchangeComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_DFNONSYMMETRICEXCHANGECOMPUTE_PARAMETERS_INT8_SLICE_COUNT = 0
   ! C name (shortened to fit Fortran's 63-char id limit): CUEST_DFNONSYMMETRICEXCHANGECOMPUTE_PARAMETERS_INT8_MODULUS_COUNT
   integer(c_int), parameter :: CUEST_DFNONSYMMETRICEXCHANGECOMPUTE_PARAM_INT8_MODULUS_COUNT = 1

   ! ---- cuestDFSymmetricDerivativeComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_DFSYMMETRICDERIVATIVECOMPUTE_PARAMETERS_MEMORY_POLICY = 0
   integer(c_int), parameter :: CUEST_DFSYMMETRICDERIVATIVECOMPUTE_PARAMETERS_FFLOAT_USAGE_MODE = 1
   integer(c_int), parameter :: CUEST_DFSYMMETRICDERIVATIVECOMPUTE_PARAMETERS_JIT_USAGE_MODE = 2

   ! ---- cuestDFSymmetricDerivativeComputeMemoryPolicy_t
   integer(c_int), parameter :: CUEST_DFSYMMETRICDERIVATIVECOMPUTE_MEMORY_POLICY_DEVICECACHE = 0
   integer(c_int), parameter :: CUEST_DFSYMMETRICDERIVATIVECOMPUTE_MEMORY_POLICY_HOSTCACHE = 1
   integer(c_int), parameter :: CUEST_DFSYMMETRICDERIVATIVECOMPUTE_MEMORY_POLICY_OVERWRITE = 2
   integer(c_int), parameter :: CUEST_DFSYMMETRICDERIVATIVECOMPUTE_MEMORY_POLICY_RECOMPUTE = 3

   ! ---- cuestXCPotentialRKSComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_XCPOTENTIALRKSCOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestXCPotentialUKSComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_XCPOTENTIALUKSCOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestXCDerivativeRKSComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_XCDERIVATIVERKSCOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestXCDerivativeUKSComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_XCDERIVATIVEUKSCOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestNonlocalXCPotentialRKSComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS_VV10_SCALE = 0
   integer(c_int), parameter :: CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS_VV10_C = 1
   integer(c_int), parameter :: CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS_VV10_B = 2
   ! C name (shortened to fit Fortran's 63-char id limit): CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAMETERS_FFLOAT_USAGE_MODE
   integer(c_int), parameter :: CUEST_NONLOCALXCPOTENTIALRKSCOMPUTE_PARAM_FFLOAT_USAGE_MODE = 3

   ! ---- cuestNonlocalXCPotentialUKSComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_NONLOCALXCPOTENTIALUKSCOMPUTE_PARAMETERS_VV10_SCALE = 0
   integer(c_int), parameter :: CUEST_NONLOCALXCPOTENTIALUKSCOMPUTE_PARAMETERS_VV10_C = 1
   integer(c_int), parameter :: CUEST_NONLOCALXCPOTENTIALUKSCOMPUTE_PARAMETERS_VV10_B = 2
   ! C name (shortened to fit Fortran's 63-char id limit): CUEST_NONLOCALXCPOTENTIALUKSCOMPUTE_PARAMETERS_FFLOAT_USAGE_MODE
   integer(c_int), parameter :: CUEST_NONLOCALXCPOTENTIALUKSCOMPUTE_PARAM_FFLOAT_USAGE_MODE = 3

   ! ---- cuestNonlocalXCDerivativeRKSComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_NONLOCALXCDERIVATIVERKSCOMPUTE_PARAMETERS_VV10_SCALE = 0
   integer(c_int), parameter :: CUEST_NONLOCALXCDERIVATIVERKSCOMPUTE_PARAMETERS_VV10_C = 1
   integer(c_int), parameter :: CUEST_NONLOCALXCDERIVATIVERKSCOMPUTE_PARAMETERS_VV10_B = 2
   ! C name (shortened to fit Fortran's 63-char id limit): CUEST_NONLOCALXCDERIVATIVERKSCOMPUTE_PARAMETERS_FFLOAT_USAGE_MODE
   integer(c_int), parameter :: CUEST_NONLOCALXCDERIVATIVERKSCOMPUTE_PARAM_FFLOAT_USAGE_MODE = 3

   ! ---- cuestNonlocalXCDerivativeUKSComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_NONLOCALXCDERIVATIVEUKSCOMPUTE_PARAMETERS_VV10_SCALE = 0
   integer(c_int), parameter :: CUEST_NONLOCALXCDERIVATIVEUKSCOMPUTE_PARAMETERS_VV10_C = 1
   integer(c_int), parameter :: CUEST_NONLOCALXCDERIVATIVEUKSCOMPUTE_PARAMETERS_VV10_B = 2
   ! C name (shortened to fit Fortran's 63-char id limit): CUEST_NONLOCALXCDERIVATIVEUKSCOMPUTE_PARAMETERS_FFLOAT_USAGE_MODE
   integer(c_int), parameter :: CUEST_NONLOCALXCDERIVATIVEUKSCOMPUTE_PARAM_FFLOAT_USAGE_MODE = 3

   ! ---- cuestXCDensityComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_XCDENSITYCOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestXCPotentialComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_XCPOTENTIALCOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestXCDerivativeComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_XCDERIVATIVECOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestXCNonsymmetricDensityComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_XCNONSYMMETRICDENSITYCOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestXCNonsymmetricDerivativeComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_XCNONSYMMETRICDERIVATIVECOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestXCGridDerivativeComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_XCGRIDDERIVATIVECOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestXCIntegrationGridComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_XCINTEGRATIONGRIDCOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestXCIntegrationWeightComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_XCINTEGRATIONWEIGHTCOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestXCIntegrationWeightComputeParametersWeightType_t
   integer(c_int), parameter :: CUEST_XCINTEGRATIONWEIGHT_PARAMETERS_WEIGHTTYPE_TOTAL = 0
   integer(c_int), parameter :: CUEST_XCINTEGRATIONWEIGHT_PARAMETERS_WEIGHTTYPE_BECKE = 1
   integer(c_int), parameter :: CUEST_XCINTEGRATIONWEIGHT_PARAMETERS_WEIGHTTYPE_QUADRATURE = 2

   ! ---- cuestXCAdvancedComputeParametersApproximation_t
   integer(c_int), parameter :: CUEST_XCADVANCED_PARAMETERS_APPROXIMATION_LDA = 0
   integer(c_int), parameter :: CUEST_XCADVANCED_PARAMETERS_APPROXIMATION_GGA = 1
   integer(c_int), parameter :: CUEST_XCADVANCED_PARAMETERS_APPROXIMATION_METAGGA = 2

   ! ---- cuestPCMIntPlanParametersAttributes_t
   integer(c_int), parameter :: CUEST_PCMINTPLAN_PARAMETERS_CUTOFF = 0
   integer(c_int), parameter :: CUEST_PCMINTPLAN_PARAMETERS_X_PREFACTOR = 1

   ! ---- cuestPCMPotentialComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_PCMPOTENTIALCOMPUTE_PARAMETERS_CONVERGENCE_THRESHOLD = 0
   integer(c_int), parameter :: CUEST_PCMPOTENTIALCOMPUTE_PARAMETERS_MAX_ITERATIONS = 1
   integer(c_int), parameter :: CUEST_PCMPOTENTIALCOMPUTE_PARAMETERS_FFLOAT_USAGE_MODE = 2
   integer(c_int), parameter :: CUEST_PCMPOTENTIALCOMPUTE_PARAMETERS_JIT_USAGE_MODE = 3

   ! ---- cuestPCMDerivativeComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_PCMDERIVATIVECOMPUTE_PARAMETERS_CONVERGENCE_THRESHOLD = 0
   integer(c_int), parameter :: CUEST_PCMDERIVATIVECOMPUTE_PARAMETERS_MAX_ITERATIONS = 1
   integer(c_int), parameter :: CUEST_PCMDERIVATIVECOMPUTE_PARAMETERS_FFLOAT_USAGE_MODE = 2
   integer(c_int), parameter :: CUEST_PCMDERIVATIVECOMPUTE_PARAMETERS_JIT_USAGE_MODE = 3

   ! ---- cuestPCMRadiiDerivativeComputeParametersAttributes_t
   ! C name (shortened to fit Fortran's 63-char id limit): CUEST_PCMRADIIDERIVATIVECOMPUTE_PARAMETERS_CONVERGENCE_THRESHOLD
   integer(c_int), parameter :: CUEST_PCMRADIIDERIVATIVECOMPUTE_PARAM_CONVERGENCE_THRESHOLD = 0
   integer(c_int), parameter :: CUEST_PCMRADIIDERIVATIVECOMPUTE_PARAMETERS_MAX_ITERATIONS = 1
   integer(c_int), parameter :: CUEST_PCMRADIIDERIVATIVECOMPUTE_PARAMETERS_FFLOAT_USAGE_MODE = 2
   integer(c_int), parameter :: CUEST_PCMRADIIDERIVATIVECOMPUTE_PARAMETERS_JIT_USAGE_MODE = 3

   ! ---- cuestMultipoleComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_MULTIPOLECOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestMultipoleDerivativeComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_MULTIPOLEDERIVATIVECOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestAngularMomentumComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_ANGULARMOMENTUMCOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestAngularMomentumDerivativeComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_ANGULARMOMENTUMDERIVATIVECOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestAngularMomentumComputeParametersComponent_t
   integer(c_int), parameter :: CUEST_ANGULARMOMENTUMCOMPUTE_PARAMETERS_COMPONENT_LX = 0
   integer(c_int), parameter :: CUEST_ANGULARMOMENTUMCOMPUTE_PARAMETERS_COMPONENT_LY = 1
   integer(c_int), parameter :: CUEST_ANGULARMOMENTUMCOMPUTE_PARAMETERS_COMPONENT_LZ = 2

   ! ---- cuestNablaComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_NABLACOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestNablaDerivativeComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_NABLADERIVATIVECOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestNablaComputeParametersComponent_t
   integer(c_int), parameter :: CUEST_NABLACOMPUTE_PARAMETERS_COMPONENT_X = 0
   integer(c_int), parameter :: CUEST_NABLACOMPUTE_PARAMETERS_COMPONENT_Y = 1
   integer(c_int), parameter :: CUEST_NABLACOMPUTE_PARAMETERS_COMPONENT_Z = 2

   ! ---- cuestPCMIntegrationGridComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_PCMINTEGRATIONGRIDCOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestPCMIntegrationWeightComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_PCMINTEGRATIONWEIGHTCOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestPCMIntegrationWeightComputeParametersWeightType_t
   integer(c_int), parameter :: CUEST_PCMINTEGRATIONWEIGHT_PARAMETERS_WEIGHTTYPE_SWITCHING = 0
   integer(c_int), parameter :: CUEST_PCMINTEGRATIONWEIGHT_PARAMETERS_WEIGHTTYPE_ZETA = 1

   ! ---- cuestDFMOIntegralsComputeParametersAttributes_t
   integer(c_int), parameter :: CUEST_DFMOINTEGRALSCOMPUTE_PARAMETERS_NONE = 0

   ! ---- cuestResultsType_t
   integer(c_int), parameter :: CUEST_PCM_RESULTS = 0

   ! ---- cuestPCMResultAttributes_t
   integer(c_int), parameter :: CUEST_PCMRESULT_PCM_DIELECTRIC_ENERGY = 0
   integer(c_int), parameter :: CUEST_PCMRESULT_CONVERGED_RESIDUAL = 1
   integer(c_int), parameter :: CUEST_PCMRESULT_NUM_ITERATIONS_TAKEN = 2
   integer(c_int), parameter :: CUEST_PCMRESULT_CONVERGED = 3

   ! ======================================================================
   !  Interoperable structs (workspace_api.h)
   ! ======================================================================
   type, bind(C) :: cuestWorkspace_t
      integer(c_intptr_t) :: hostBuffer = 0_c_intptr_t
      integer(c_size_t)   :: hostBufferSizeInBytes = 0_c_size_t
      integer(c_intptr_t) :: deviceBuffer = 0_c_intptr_t
      integer(c_size_t)   :: deviceBufferSizeInBytes = 0_c_size_t
   end type cuestWorkspace_t

   type, bind(C) :: cuestWorkspaceDescriptor_t
      integer(c_size_t)   :: hostBufferSizeInBytes = 0_c_size_t
      integer(c_size_t)   :: deviceBufferSizeInBytes = 0_c_size_t
   end type cuestWorkspaceDescriptor_t

   ! ======================================================================
   !  C function interfaces
   ! ======================================================================
   interface

      ! ---------------------------------------------------------------
      !  cuest/context/context_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestCreate(parameters, handle) &
         bind(C, name="cuestCreate")
         import
         type(c_ptr), value :: parameters
         type(c_ptr), intent(out) :: handle
      end function cuestCreate

      integer(c_int) function cuestDestroy(handle) &
         bind(C, name="cuestDestroy")
         import
         type(c_ptr), value :: handle
      end function cuestDestroy

      integer(c_int) function cuestSetMathMode(handle, mode) &
         bind(C, name="cuestSetMathMode")
         import
         type(c_ptr), value :: handle
         integer(c_int), value :: mode
      end function cuestSetMathMode

      integer(c_int) function cuestGetMathMode(handle, mode) &
         bind(C, name="cuestGetMathMode")
         import
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: mode
      end function cuestGetMathMode

      integer(c_int) function cuestSetComputeCapabilityTarget(handle, targetMajorVersion, targetMinorVersion) &
         bind(C, name="cuestSetComputeCapabilityTarget")
         import
         type(c_ptr), value :: handle
         integer(c_int32_t), value :: targetMajorVersion
         integer(c_int32_t), value :: targetMinorVersion
      end function cuestSetComputeCapabilityTarget

      integer(c_int) function cuestGetComputeCapabilityTarget(handle, targetMajorVersion, targetMinorVersion) &
         bind(C, name="cuestGetComputeCapabilityTarget")
         import
         type(c_ptr), value :: handle
         integer(c_int32_t), intent(out) :: targetMajorVersion
         integer(c_int32_t), intent(out) :: targetMinorVersion
      end function cuestGetComputeCapabilityTarget

      integer(c_int) function cuestGetMajorVersion(handle, majorVersion) &
         bind(C, name="cuestGetMajorVersion")
         import
         type(c_ptr), value :: handle
         integer(c_int32_t), intent(out) :: majorVersion
      end function cuestGetMajorVersion

      integer(c_int) function cuestGetMinorVersion(handle, minorVersion) &
         bind(C, name="cuestGetMinorVersion")
         import
         type(c_ptr), value :: handle
         integer(c_int32_t), intent(out) :: minorVersion
      end function cuestGetMinorVersion

      integer(c_int) function cuestGetPatchVersion(handle, patchVersion) &
         bind(C, name="cuestGetPatchVersion")
         import
         type(c_ptr), value :: handle
         integer(c_int32_t), intent(out) :: patchVersion
      end function cuestGetPatchVersion

      ! ---------------------------------------------------------------
      !  cuest/basis/ao_shell_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestAOShellCreate( &
         handle, isPure, L, numPrimitive, exponents, coefficients, parameters, outShell) &
         bind(C, name="cuestAOShellCreate")
         import
         type(c_ptr), value :: handle
         integer(c_int32_t), value :: isPure
         integer(c_int64_t), value :: L
         integer(c_int64_t), value :: numPrimitive
         type(c_ptr), value :: exponents
         type(c_ptr), value :: coefficients
         type(c_ptr), value :: parameters
         type(c_ptr), intent(out) :: outShell
      end function cuestAOShellCreate

      integer(c_int) function cuestAOShellDestroy(shell) &
         bind(C, name="cuestAOShellDestroy")
         import
         type(c_ptr), value :: shell
      end function cuestAOShellDestroy

      ! ---------------------------------------------------------------
      !  cuest/basis/ao_basis_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestAOBasisCreate( &
         handle, numAtoms, numShellsPerAtom, shells, parameters, persistentWorkspace, temporaryWorkspace, outBasis) &
         bind(C, name="cuestAOBasisCreate")
         import
         type(c_ptr), value :: handle
         integer(c_int64_t), value :: numAtoms
         integer(c_int64_t), dimension(*), intent(in) :: numShellsPerAtom
         type(c_ptr), dimension(*), intent(in) :: shells
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: persistentWorkspace
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), intent(out) :: outBasis
      end function cuestAOBasisCreate

      integer(c_int) function cuestAOBasisCreateWorkspaceQuery( &
         handle, numAtoms, numShellsPerAtom, shells, parameters, persistentWorkspaceDescriptor, &
         temporaryWorkspaceDescriptor, outBasis) &
         bind(C, name="cuestAOBasisCreateWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         integer(c_int64_t), value :: numAtoms
         integer(c_int64_t), dimension(*), intent(in) :: numShellsPerAtom
         type(c_ptr), dimension(*), intent(in) :: shells
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: persistentWorkspaceDescriptor
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), intent(out) :: outBasis
      end function cuestAOBasisCreateWorkspaceQuery

      integer(c_int) function cuestAOBasisDestroy(basis) &
         bind(C, name="cuestAOBasisDestroy")
         import
         type(c_ptr), value :: basis
      end function cuestAOBasisDestroy

      ! ---------------------------------------------------------------
      !  cuest/basis/ao_pair_list_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestAOPairListCreate( &
         handle, basis, numAtoms, xyz, thresholdPQ, parameters, persistentWorkspace, temporaryWorkspace, &
         outPairList) &
         bind(C, name="cuestAOPairListCreate")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: basis
         integer(c_int64_t), value :: numAtoms
         type(c_ptr), value :: xyz
         real(c_double), value :: thresholdPQ
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: persistentWorkspace
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), intent(out) :: outPairList
      end function cuestAOPairListCreate

      integer(c_int) function cuestAOPairListCreateWorkspaceQuery( &
         handle, basis, numAtoms, xyz, thresholdPQ, parameters, persistentWorkspaceDescriptor, &
         temporaryWorkspaceDescriptor, outPairList) &
         bind(C, name="cuestAOPairListCreateWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: basis
         integer(c_int64_t), value :: numAtoms
         type(c_ptr), value :: xyz
         real(c_double), value :: thresholdPQ
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: persistentWorkspaceDescriptor
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), intent(out) :: outPairList
      end function cuestAOPairListCreateWorkspaceQuery

      integer(c_int) function cuestAOPairListDestroy(pairList) &
         bind(C, name="cuestAOPairListDestroy")
         import
         type(c_ptr), value :: pairList
      end function cuestAOPairListDestroy

      ! ---------------------------------------------------------------
      !  cuest/ecp/ecp_integral_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestECPCompute( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspace, outECPMatrix) &
         bind(C, name="cuestECPCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), value :: outECPMatrix
      end function cuestECPCompute

      integer(c_int) function cuestECPComputeWorkspaceQuery( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspaceDescriptor, outECPMatrix) &
         bind(C, name="cuestECPComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), value :: outECPMatrix
      end function cuestECPComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/ecp/ecp_integral_derivative_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestECPDerivativeCompute( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspace, densityMatrix, outGradient) &
         bind(C, name="cuestECPDerivativeCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: outGradient
      end function cuestECPDerivativeCompute

      integer(c_int) function cuestECPDerivativeComputeWorkspaceQuery( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspaceDescriptor, densityMatrix, outGradient) &
         bind(C, name="cuestECPDerivativeComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: outGradient
      end function cuestECPDerivativeComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/ecp_basis/ecp_atom_api.h
      ! ---------------------------------------------------------------
   integer(c_int) function cuestECPAtomCreate(handle, numElectrons, numShells, shells, topShell, parameters, outECPAtom) &
         bind(C, name="cuestECPAtomCreate")
         import
         type(c_ptr), value :: handle
         integer(c_int64_t), value :: numElectrons
         integer(c_int64_t), value :: numShells
         type(c_ptr), dimension(*), intent(in) :: shells
         type(c_ptr), value :: topShell
         type(c_ptr), value :: parameters
         type(c_ptr), intent(out) :: outECPAtom
      end function cuestECPAtomCreate

      integer(c_int) function cuestECPAtomDestroy(atom) &
         bind(C, name="cuestECPAtomDestroy")
         import
         type(c_ptr), value :: atom
      end function cuestECPAtomDestroy

      ! ---------------------------------------------------------------
      !  cuest/ecp_basis/ecp_shell_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestECPShellCreate( &
         handle, L, numPrimitive, radialPowers, coefficients, exponents, parameters, outECPShell) &
         bind(C, name="cuestECPShellCreate")
         import
         type(c_ptr), value :: handle
         integer(c_int64_t), value :: L
         integer(c_int64_t), value :: numPrimitive
         integer(c_int64_t), dimension(*), intent(in) :: radialPowers
         type(c_ptr), value :: coefficients
         type(c_ptr), value :: exponents
         type(c_ptr), value :: parameters
         type(c_ptr), intent(out) :: outECPShell
      end function cuestECPShellCreate

      integer(c_int) function cuestECPShellDestroy(shell) &
         bind(C, name="cuestECPShellDestroy")
         import
         type(c_ptr), value :: shell
      end function cuestECPShellDestroy

      ! ---------------------------------------------------------------
      !  cuest/integral_plans/ecp_integral_plan_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestECPIntPlanCreate( &
         handle, basis, xyzCPU, numECPAtoms, activeIndices, activeAtoms, parameters, persistentWorkspace, &
         temporaryWorkspace, outPlan) &
         bind(C, name="cuestECPIntPlanCreate")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: basis
         type(c_ptr), value :: xyzCPU
         integer(c_int64_t), value :: numECPAtoms
         integer(c_int64_t), dimension(*), intent(in) :: activeIndices
         type(c_ptr), dimension(*), intent(in) :: activeAtoms
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: persistentWorkspace
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), intent(out) :: outPlan
      end function cuestECPIntPlanCreate

      integer(c_int) function cuestECPIntPlanCreateWorkspaceQuery( &
         handle, basis, xyzCPU, numECPAtoms, activeIndices, activeAtoms, parameters, persistentWorkspaceDescriptor, &
         temporaryWorkspaceDescriptor, outPlan) &
         bind(C, name="cuestECPIntPlanCreateWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: basis
         type(c_ptr), value :: xyzCPU
         integer(c_int64_t), value :: numECPAtoms
         integer(c_int64_t), dimension(*), intent(in) :: activeIndices
         type(c_ptr), dimension(*), intent(in) :: activeAtoms
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: persistentWorkspaceDescriptor
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), intent(out) :: outPlan
      end function cuestECPIntPlanCreateWorkspaceQuery

      integer(c_int) function cuestECPIntPlanDestroy(plan) &
         bind(C, name="cuestECPIntPlanDestroy")
         import
         type(c_ptr), value :: plan
      end function cuestECPIntPlanDestroy

      ! ---------------------------------------------------------------
      !  cuest/integral_plans/one_electron_integral_plan_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestOEIntPlanCreate( &
         handle, basis, pairList, parameters, persistentWorkspace, temporaryWorkspace, outPlan) &
         bind(C, name="cuestOEIntPlanCreate")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: basis
         type(c_ptr), value :: pairList
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: persistentWorkspace
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), intent(out) :: outPlan
      end function cuestOEIntPlanCreate

      integer(c_int) function cuestOEIntPlanCreateWorkspaceQuery( &
         handle, basis, pairList, parameters, persistentWorkspaceDescriptor, temporaryWorkspaceDescriptor, outPlan) &
         bind(C, name="cuestOEIntPlanCreateWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: basis
         type(c_ptr), value :: pairList
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: persistentWorkspaceDescriptor
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), intent(out) :: outPlan
      end function cuestOEIntPlanCreateWorkspaceQuery

      integer(c_int) function cuestOEIntPlanDestroy(plan) &
         bind(C, name="cuestOEIntPlanDestroy")
         import
         type(c_ptr), value :: plan
      end function cuestOEIntPlanDestroy

      ! ---------------------------------------------------------------
      !  cuest/integral_plans/density_fitted_integral_plan_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestDFIntPlanCreate( &
         handle, primaryBasis, auxiliaryBasis, pairList, parameters, persistentWorkspace, temporaryWorkspace, &
         outPlan) &
         bind(C, name="cuestDFIntPlanCreate")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: primaryBasis
         type(c_ptr), value :: auxiliaryBasis
         type(c_ptr), value :: pairList
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: persistentWorkspace
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), intent(out) :: outPlan
      end function cuestDFIntPlanCreate

      integer(c_int) function cuestDFIntPlanCreateWorkspaceQuery( &
         handle, primaryBasis, auxiliaryBasis, pairList, parameters, persistentWorkspaceDescriptor, &
         temporaryWorkspaceDescriptor, outPlan) &
         bind(C, name="cuestDFIntPlanCreateWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: primaryBasis
         type(c_ptr), value :: auxiliaryBasis
         type(c_ptr), value :: pairList
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: persistentWorkspaceDescriptor
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), intent(out) :: outPlan
      end function cuestDFIntPlanCreateWorkspaceQuery

      integer(c_int) function cuestDFIntPlanDestroy(plan) &
         bind(C, name="cuestDFIntPlanDestroy")
         import
         type(c_ptr), value :: plan
      end function cuestDFIntPlanDestroy

      ! ---------------------------------------------------------------
      !  cuest/one_electron_integrals/overlap_integral_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestOverlapCompute(handle, plan, parameters, temporaryWorkspace, outSMatrix) &
         bind(C, name="cuestOverlapCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), value :: outSMatrix
      end function cuestOverlapCompute

      integer(c_int) function cuestOverlapComputeWorkspaceQuery( &
         handle, plan, parameters, temporaryWorkspaceDescriptor, outSMatrix) &
         bind(C, name="cuestOverlapComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), value :: outSMatrix
      end function cuestOverlapComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/one_electron_integrals/overlap_integral_derivative_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestOverlapDerivativeCompute( &
         handle, plan, parameters, temporaryWorkspace, densityMatrix, outGradient) &
         bind(C, name="cuestOverlapDerivativeCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: outGradient
      end function cuestOverlapDerivativeCompute

      integer(c_int) function cuestOverlapDerivativeComputeWorkspaceQuery( &
         handle, plan, parameters, temporaryWorkspaceDescriptor, densityMatrix, outGradient) &
         bind(C, name="cuestOverlapDerivativeComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: outGradient
      end function cuestOverlapDerivativeComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/one_electron_integrals/kinetic_integral_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestKineticCompute(handle, plan, parameters, temporaryWorkspace, outTMatrix) &
         bind(C, name="cuestKineticCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), value :: outTMatrix
      end function cuestKineticCompute

      integer(c_int) function cuestKineticComputeWorkspaceQuery( &
         handle, plan, parameters, temporaryWorkspaceDescriptor, outTMatrix) &
         bind(C, name="cuestKineticComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), value :: outTMatrix
      end function cuestKineticComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/one_electron_integrals/kinetic_integral_derivative_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestKineticDerivativeCompute( &
         handle, plan, parameters, temporaryWorkspace, densityMatrix, outGradient) &
         bind(C, name="cuestKineticDerivativeCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: outGradient
      end function cuestKineticDerivativeCompute

      integer(c_int) function cuestKineticDerivativeComputeWorkspaceQuery( &
         handle, plan, parameters, temporaryWorkspaceDescriptor, densityMatrix, outGradient) &
         bind(C, name="cuestKineticDerivativeComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: outGradient
      end function cuestKineticDerivativeComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/one_electron_integrals/multipole_integral_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestMultipoleCompute( &
         handle, plan, parameters, temporaryWorkspace, multipoleOrder, origin, outMatrix) &
         bind(C, name="cuestMultipoleCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int32_t), dimension(*), intent(in) :: multipoleOrder
         type(c_ptr), value :: origin
         type(c_ptr), value :: outMatrix
      end function cuestMultipoleCompute

      integer(c_int) function cuestMultipoleComputeWorkspaceQuery( &
         handle, plan, parameters, temporaryWorkspaceDescriptor, multipoleOrder, origin, outMatrix) &
         bind(C, name="cuestMultipoleComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int32_t), dimension(*), intent(in) :: multipoleOrder
         type(c_ptr), value :: origin
         type(c_ptr), value :: outMatrix
      end function cuestMultipoleComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/one_electron_integrals/multipole_integral_derivative_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestMultipoleDerivativeCompute( &
         handle, plan, parameters, temporaryWorkspace, multipoleOrder, origin, densityMatrix, outGradient) &
         bind(C, name="cuestMultipoleDerivativeCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int32_t), dimension(*), intent(in) :: multipoleOrder
         type(c_ptr), value :: origin
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: outGradient
      end function cuestMultipoleDerivativeCompute

      integer(c_int) function cuestMultipoleDerivativeComputeWorkspaceQuery( &
         handle, plan, parameters, temporaryWorkspaceDescriptor, multipoleOrder, origin, densityMatrix, outGradient) &
         bind(C, name="cuestMultipoleDerivativeComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int32_t), dimension(*), intent(in) :: multipoleOrder
         type(c_ptr), value :: origin
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: outGradient
      end function cuestMultipoleDerivativeComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/one_electron_integrals/angular_momentum_integral_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestAngularMomentumCompute( &
         handle, plan, parameters, temporaryWorkspace, component, origin, outMatrix) &
         bind(C, name="cuestAngularMomentumCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int), value :: component
         type(c_ptr), value :: origin
         type(c_ptr), value :: outMatrix
      end function cuestAngularMomentumCompute

      integer(c_int) function cuestAngularMomentumComputeWorkspaceQuery( &
         handle, plan, parameters, temporaryWorkspaceDescriptor, component, origin, outMatrix) &
         bind(C, name="cuestAngularMomentumComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int), value :: component
         type(c_ptr), value :: origin
         type(c_ptr), value :: outMatrix
      end function cuestAngularMomentumComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/one_electron_integrals/angular_momentum_integral_derivative_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestAngularMomentumDerivativeCompute( &
         handle, plan, parameters, temporaryWorkspace, component, origin, densityMatrix, outGradient) &
         bind(C, name="cuestAngularMomentumDerivativeCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int), value :: component
         type(c_ptr), value :: origin
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: outGradient
      end function cuestAngularMomentumDerivativeCompute

      integer(c_int) function cuestAngularMomentumDerivativeComputeWorkspaceQuery( &
         handle, plan, parameters, temporaryWorkspaceDescriptor, component, origin, densityMatrix, outGradient) &
         bind(C, name="cuestAngularMomentumDerivativeComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int), value :: component
         type(c_ptr), value :: origin
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: outGradient
      end function cuestAngularMomentumDerivativeComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/one_electron_integrals/nabla_integral_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestNablaCompute(handle, plan, parameters, temporaryWorkspace, component, outMatrix) &
         bind(C, name="cuestNablaCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int), value :: component
         type(c_ptr), value :: outMatrix
      end function cuestNablaCompute

      integer(c_int) function cuestNablaComputeWorkspaceQuery( &
         handle, plan, parameters, temporaryWorkspaceDescriptor, component, outMatrix) &
         bind(C, name="cuestNablaComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int), value :: component
         type(c_ptr), value :: outMatrix
      end function cuestNablaComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/one_electron_integrals/nabla_integral_derivative_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestNablaDerivativeCompute( &
         handle, plan, parameters, temporaryWorkspace, component, densityMatrix, outGradient) &
         bind(C, name="cuestNablaDerivativeCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int), value :: component
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: outGradient
      end function cuestNablaDerivativeCompute

      integer(c_int) function cuestNablaDerivativeComputeWorkspaceQuery( &
         handle, plan, parameters, temporaryWorkspaceDescriptor, component, densityMatrix, outGradient) &
         bind(C, name="cuestNablaDerivativeComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int), value :: component
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: outGradient
      end function cuestNablaDerivativeComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/one_electron_integrals/potential_integral_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestPotentialCompute( &
         handle, plan, parameters, temporaryWorkspace, numCharges, xyz, q, outVMatrix) &
         bind(C, name="cuestPotentialCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int64_t), value :: numCharges
         type(c_ptr), value :: xyz
         type(c_ptr), value :: q
         type(c_ptr), value :: outVMatrix
      end function cuestPotentialCompute

      integer(c_int) function cuestPotentialComputeWorkspaceQuery( &
         handle, plan, parameters, temporaryWorkspaceDescriptor, numCharges, xyz, q, outVMatrix) &
         bind(C, name="cuestPotentialComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int64_t), value :: numCharges
         type(c_ptr), value :: xyz
         type(c_ptr), value :: q
         type(c_ptr), value :: outVMatrix
      end function cuestPotentialComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/one_electron_integrals/potential_integral_derivative_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestPotentialDerivativeCompute( &
         handle, plan, parameters, temporaryWorkspace, numCharges, xyz, q, densityMatrix, outBasisGradient, &
         outChargeGradient) &
         bind(C, name="cuestPotentialDerivativeCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int64_t), value :: numCharges
         type(c_ptr), value :: xyz
         type(c_ptr), value :: q
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: outBasisGradient
         type(c_ptr), value :: outChargeGradient
      end function cuestPotentialDerivativeCompute

      integer(c_int) function cuestPotentialDerivativeComputeWorkspaceQuery( &
         handle, plan, parameters, temporaryWorkspaceDescriptor, numCharges, xyz, q, densityMatrix, &
         outBasisGradient, outChargeGradient) &
         bind(C, name="cuestPotentialDerivativeComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int64_t), value :: numCharges
         type(c_ptr), value :: xyz
         type(c_ptr), value :: q
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: outBasisGradient
         type(c_ptr), value :: outChargeGradient
      end function cuestPotentialDerivativeComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/density_fitted_integrals/density_fitted_coulomb_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestDFCoulombCompute( &
         handle, plan, parameters, temporaryWorkspace, densityMatrix, outCoulombMatrix) &
         bind(C, name="cuestDFCoulombCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: outCoulombMatrix
      end function cuestDFCoulombCompute

      integer(c_int) function cuestDFCoulombComputeWorkspaceQuery( &
         handle, plan, parameters, temporaryWorkspaceDescriptor, densityMatrix, outCoulombMatrix) &
         bind(C, name="cuestDFCoulombComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: outCoulombMatrix
      end function cuestDFCoulombComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/density_fitted_integrals/density_fitted_exchange_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestDFSymmetricExchangeCompute( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspace, numOccupied, coefficientMatrix, &
         outExchangeMatrix) &
         bind(C, name="cuestDFSymmetricExchangeCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: coefficientMatrix
         type(c_ptr), value :: outExchangeMatrix
      end function cuestDFSymmetricExchangeCompute

      integer(c_int) function cuestDFSymmetricExchangeComputeWorkspaceQuery( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspaceDescriptor, numOccupied, &
         coefficientMatrix, outExchangeMatrix) &
         bind(C, name="cuestDFSymmetricExchangeComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: coefficientMatrix
         type(c_ptr), value :: outExchangeMatrix
      end function cuestDFSymmetricExchangeComputeWorkspaceQuery

      integer(c_int) function cuestDFNonsymmetricExchangeCompute( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspace, numCoefficientMatrices, numOccupied, &
         leftCoefficientMatrix, rightCoefficientMatrices, outExchangeMatrices) &
         bind(C, name="cuestDFNonsymmetricExchangeCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int64_t), value :: numCoefficientMatrices
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: leftCoefficientMatrix
         type(c_ptr), value :: rightCoefficientMatrices
         type(c_ptr), value :: outExchangeMatrices
      end function cuestDFNonsymmetricExchangeCompute

      integer(c_int) function cuestDFNonsymmetricExchangeComputeWorkspaceQuery( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspaceDescriptor, numCoefficientMatrices, &
         numOccupied, leftCoefficientMatrix, rightCoefficientMatrices, outExchangeMatrices) &
         bind(C, name="cuestDFNonsymmetricExchangeComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int64_t), value :: numCoefficientMatrices
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: leftCoefficientMatrix
         type(c_ptr), value :: rightCoefficientMatrices
         type(c_ptr), value :: outExchangeMatrices
      end function cuestDFNonsymmetricExchangeComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/density_fitted_integrals/density_fitted_derivative_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestDFSymmetricDerivativeCompute( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspace, densityScale, densityMatrix, &
         coefficientScale, numCoefficientMatrices, numOccupied, coefficientMatrices, outGradient) &
         bind(C, name="cuestDFSymmetricDerivativeCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspace_t) :: temporaryWorkspace
         real(c_double), value :: densityScale
         type(c_ptr), value :: densityMatrix
         real(c_double), value :: coefficientScale
         integer(c_int64_t), value :: numCoefficientMatrices
         integer(c_int64_t), dimension(*), intent(in) :: numOccupied
         type(c_ptr), value :: coefficientMatrices
         type(c_ptr), value :: outGradient
      end function cuestDFSymmetricDerivativeCompute

      integer(c_int) function cuestDFSymmetricDerivativeComputeWorkspaceQuery( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspaceDescriptor, densityScale, densityMatrix, &
         coefficientScale, numCoefficientMatrices, numOccupied, coefficientMatrices, outGradient) &
         bind(C, name="cuestDFSymmetricDerivativeComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         real(c_double), value :: densityScale
         type(c_ptr), value :: densityMatrix
         real(c_double), value :: coefficientScale
         integer(c_int64_t), value :: numCoefficientMatrices
         integer(c_int64_t), dimension(*), intent(in) :: numOccupied
         type(c_ptr), value :: coefficientMatrices
         type(c_ptr), value :: outGradient
      end function cuestDFSymmetricDerivativeComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/density_fitted_integrals/density_fitted_mo_integrals_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestDFMOIntegralsCompute( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspace, numCoefficientMatrices, numLeftOrbitals, &
         numRightOrbitals, leftCoefficientMatrices, rightCoefficientMatrices, outTensors) &
         bind(C, name="cuestDFMOIntegralsCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int64_t), value :: numCoefficientMatrices
         integer(c_int64_t), dimension(*), intent(in) :: numLeftOrbitals
         integer(c_int64_t), dimension(*), intent(in) :: numRightOrbitals
         type(c_ptr), value :: leftCoefficientMatrices
         type(c_ptr), value :: rightCoefficientMatrices
         type(c_ptr), value :: outTensors
      end function cuestDFMOIntegralsCompute

      integer(c_int) function cuestDFMOIntegralsComputeWorkspaceQuery( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspaceDescriptor, numCoefficientMatrices, &
         numLeftOrbitals, numRightOrbitals, leftCoefficientMatrices, rightCoefficientMatrices, outTensors) &
         bind(C, name="cuestDFMOIntegralsComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int64_t), value :: numCoefficientMatrices
         integer(c_int64_t), dimension(*), intent(in) :: numLeftOrbitals
         integer(c_int64_t), dimension(*), intent(in) :: numRightOrbitals
         type(c_ptr), value :: leftCoefficientMatrices
         type(c_ptr), value :: rightCoefficientMatrices
         type(c_ptr), value :: outTensors
      end function cuestDFMOIntegralsComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/grid/atom_grid_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestAtomGridCreate( &
         handle, numRadialPoints, radialNodes, radialWeights, numAngularPoints, parameters, outAtomGrid) &
         bind(C, name="cuestAtomGridCreate")
         import
         type(c_ptr), value :: handle
         integer(c_int64_t), value :: numRadialPoints
         type(c_ptr), value :: radialNodes
         type(c_ptr), value :: radialWeights
         integer(c_int64_t), dimension(*), intent(in) :: numAngularPoints
         type(c_ptr), value :: parameters
         type(c_ptr), intent(out) :: outAtomGrid
      end function cuestAtomGridCreate

      integer(c_int) function cuestAtomGridDestroy(atomGrid) &
         bind(C, name="cuestAtomGridDestroy")
         import
         type(c_ptr), value :: atomGrid
      end function cuestAtomGridDestroy

      ! ---------------------------------------------------------------
      !  cuest/grid/molecular_grid_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestMolecularGridCreate( &
         handle, numAtoms, atomGrid, xyz, parameters, persistentWorkspace, temporaryWorkspace, outGrid) &
         bind(C, name="cuestMolecularGridCreate")
         import
         type(c_ptr), value :: handle
         integer(c_int64_t), value :: numAtoms
         type(c_ptr), dimension(*), intent(in) :: atomGrid
         type(c_ptr), value :: xyz
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: persistentWorkspace
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), intent(out) :: outGrid
      end function cuestMolecularGridCreate

      integer(c_int) function cuestMolecularGridCreateWorkspaceQuery( &
         handle, numAtoms, atomGrid, xyz, parameters, persistentWorkspaceDescriptor, temporaryWorkspaceDescriptor, &
         outGrid) &
         bind(C, name="cuestMolecularGridCreateWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         integer(c_int64_t), value :: numAtoms
         type(c_ptr), dimension(*), intent(in) :: atomGrid
         type(c_ptr), value :: xyz
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: persistentWorkspaceDescriptor
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), intent(out) :: outGrid
      end function cuestMolecularGridCreateWorkspaceQuery

      integer(c_int) function cuestMolecularGridDestroy(grid) &
         bind(C, name="cuestMolecularGridDestroy")
         import
         type(c_ptr), value :: grid
      end function cuestMolecularGridDestroy

      ! ---------------------------------------------------------------
      !  cuest/integral_plans/xc_integral_plan_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestXCIntPlanCreate( &
         handle, basis, grid, functional, parameters, persistentWorkspace, temporaryWorkspace, outPlan) &
         bind(C, name="cuestXCIntPlanCreate")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: basis
         type(c_ptr), value :: grid
         integer(c_int), value :: functional
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: persistentWorkspace
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), intent(out) :: outPlan
      end function cuestXCIntPlanCreate

      integer(c_int) function cuestXCIntPlanCreateWorkspaceQuery( &
         handle, basis, grid, functional, parameters, persistentWorkspaceDescriptor, temporaryWorkspaceDescriptor, &
         outPlan) &
         bind(C, name="cuestXCIntPlanCreateWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: basis
         type(c_ptr), value :: grid
         integer(c_int), value :: functional
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: persistentWorkspaceDescriptor
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), intent(out) :: outPlan
      end function cuestXCIntPlanCreateWorkspaceQuery

      integer(c_int) function cuestXCIntPlanDestroy(plan) &
         bind(C, name="cuestXCIntPlanDestroy")
         import
         type(c_ptr), value :: plan
      end function cuestXCIntPlanDestroy

      ! ---------------------------------------------------------------
      !  cuest/exchange_correlation/xc_integration_grid_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestXCIntegrationGridCompute(handle, plan, parameters, temporaryWorkspace, outGridPoints) &
         bind(C, name="cuestXCIntegrationGridCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), value :: outGridPoints
      end function cuestXCIntegrationGridCompute

      integer(c_int) function cuestXCIntegrationGridComputeWorkspaceQuery( &
         handle, plan, parameters, temporaryWorkspaceDescriptor, outGridPoints) &
         bind(C, name="cuestXCIntegrationGridComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), value :: outGridPoints
      end function cuestXCIntegrationGridComputeWorkspaceQuery

      integer(c_int) function cuestXCIntegrationWeightCompute( &
         handle, plan, weightType, parameters, temporaryWorkspace, outGridWeights) &
         bind(C, name="cuestXCIntegrationWeightCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         integer(c_int), value :: weightType
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), value :: outGridWeights
      end function cuestXCIntegrationWeightCompute

      integer(c_int) function cuestXCIntegrationWeightComputeWorkspaceQuery( &
         handle, plan, weightType, parameters, temporaryWorkspaceDescriptor, outGridWeights) &
         bind(C, name="cuestXCIntegrationWeightComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         integer(c_int), value :: weightType
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), value :: outGridWeights
      end function cuestXCIntegrationWeightComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/exchange_correlation/xc_potential_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestXCPotentialRKSCompute( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspace, numOccupied, coefficientMatrix, &
         outXCEnergy, outXCPotentialMatrix) &
         bind(C, name="cuestXCPotentialRKSCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: coefficientMatrix
         type(c_ptr), value :: outXCEnergy
         type(c_ptr), value :: outXCPotentialMatrix
      end function cuestXCPotentialRKSCompute

      integer(c_int) function cuestXCPotentialRKSComputeWorkspaceQuery( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspaceDescriptor, numOccupied, &
         coefficientMatrix, outXCEnergy, outXCPotentialMatrix) &
         bind(C, name="cuestXCPotentialRKSComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: coefficientMatrix
         type(c_ptr), value :: outXCEnergy
         type(c_ptr), value :: outXCPotentialMatrix
      end function cuestXCPotentialRKSComputeWorkspaceQuery

      integer(c_int) function cuestXCPotentialUKSCompute( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspace, numOccupiedAlpha, numOccupiedBeta, &
         coefficientMatrixAlpha, coefficientMatrixBeta, outXCEnergy, outXCPotentialMatrixAlpha, &
         outXCPotentialMatrixBeta) &
         bind(C, name="cuestXCPotentialUKSCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int64_t), value :: numOccupiedAlpha
         integer(c_int64_t), value :: numOccupiedBeta
         type(c_ptr), value :: coefficientMatrixAlpha
         type(c_ptr), value :: coefficientMatrixBeta
         type(c_ptr), value :: outXCEnergy
         type(c_ptr), value :: outXCPotentialMatrixAlpha
         type(c_ptr), value :: outXCPotentialMatrixBeta
      end function cuestXCPotentialUKSCompute

      integer(c_int) function cuestXCPotentialUKSComputeWorkspaceQuery( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspaceDescriptor, numOccupiedAlpha, &
         numOccupiedBeta, coefficientMatrixAlpha, coefficientMatrixBeta, outXCEnergy, outXCPotentialMatrixAlpha, &
         outXCPotentialMatrixBeta) &
         bind(C, name="cuestXCPotentialUKSComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int64_t), value :: numOccupiedAlpha
         integer(c_int64_t), value :: numOccupiedBeta
         type(c_ptr), value :: coefficientMatrixAlpha
         type(c_ptr), value :: coefficientMatrixBeta
         type(c_ptr), value :: outXCEnergy
         type(c_ptr), value :: outXCPotentialMatrixAlpha
         type(c_ptr), value :: outXCPotentialMatrixBeta
      end function cuestXCPotentialUKSComputeWorkspaceQuery

      integer(c_int) function cuestXCDensityCompute( &
         handle, plan, approximation, parameters, variableBufferSize, temporaryWorkspace, numOccupied, &
         coefficientMatrix, outGridDensity) &
         bind(C, name="cuestXCDensityCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         integer(c_int), value :: approximation
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: coefficientMatrix
         type(c_ptr), value :: outGridDensity
      end function cuestXCDensityCompute

      integer(c_int) function cuestXCDensityComputeWorkspaceQuery( &
         handle, plan, approximation, parameters, variableBufferSize, temporaryWorkspaceDescriptor, numOccupied, &
         coefficientMatrix, outGridDensity) &
         bind(C, name="cuestXCDensityComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         integer(c_int), value :: approximation
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: coefficientMatrix
         type(c_ptr), value :: outGridDensity
      end function cuestXCDensityComputeWorkspaceQuery

      integer(c_int) function cuestXCPotentialCompute( &
         handle, plan, approximation, parameters, variableBufferSize, temporaryWorkspace, gridXCPotential, &
         outXCPotentialMatrix) &
         bind(C, name="cuestXCPotentialCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         integer(c_int), value :: approximation
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), value :: gridXCPotential
         type(c_ptr), value :: outXCPotentialMatrix
      end function cuestXCPotentialCompute

      integer(c_int) function cuestXCPotentialComputeWorkspaceQuery( &
         handle, plan, approximation, parameters, variableBufferSize, temporaryWorkspaceDescriptor, &
         gridXCPotential, outXCPotentialMatrix) &
         bind(C, name="cuestXCPotentialComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         integer(c_int), value :: approximation
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), value :: gridXCPotential
         type(c_ptr), value :: outXCPotentialMatrix
      end function cuestXCPotentialComputeWorkspaceQuery

      integer(c_int) function cuestXCNonsymmetricDensityCompute( &
         handle, plan, approximation, parameters, variableBufferSize, temporaryWorkspace, numOccupied, &
         leftCoefficientMatrix, rightCoefficientMatrix, outGridDensity) &
         bind(C, name="cuestXCNonsymmetricDensityCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         integer(c_int), value :: approximation
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: leftCoefficientMatrix
         type(c_ptr), value :: rightCoefficientMatrix
         type(c_ptr), value :: outGridDensity
      end function cuestXCNonsymmetricDensityCompute

      integer(c_int) function cuestXCNonsymmetricDensityComputeWorkspaceQuery( &
         handle, plan, approximation, parameters, variableBufferSize, temporaryWorkspaceDescriptor, numOccupied, &
         leftCoefficientMatrix, rightCoefficientMatrix, outGridDensity) &
         bind(C, name="cuestXCNonsymmetricDensityComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         integer(c_int), value :: approximation
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: leftCoefficientMatrix
         type(c_ptr), value :: rightCoefficientMatrix
         type(c_ptr), value :: outGridDensity
      end function cuestXCNonsymmetricDensityComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/exchange_correlation/xc_derivative_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestXCDerivativeRKSCompute( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspace, numOccupied, coefficientMatrix, &
         outGradient) &
         bind(C, name="cuestXCDerivativeRKSCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: coefficientMatrix
         type(c_ptr), value :: outGradient
      end function cuestXCDerivativeRKSCompute

      integer(c_int) function cuestXCDerivativeRKSComputeWorkspaceQuery( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspaceDescriptor, numOccupied, &
         coefficientMatrix, outGradient) &
         bind(C, name="cuestXCDerivativeRKSComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: coefficientMatrix
         type(c_ptr), value :: outGradient
      end function cuestXCDerivativeRKSComputeWorkspaceQuery

      integer(c_int) function cuestXCDerivativeUKSCompute( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspace, numOccupiedAlpha, numOccupiedBeta, &
         coefficientMatrixAlpha, coefficientMatrixBeta, outGradient) &
         bind(C, name="cuestXCDerivativeUKSCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int64_t), value :: numOccupiedAlpha
         integer(c_int64_t), value :: numOccupiedBeta
         type(c_ptr), value :: coefficientMatrixAlpha
         type(c_ptr), value :: coefficientMatrixBeta
         type(c_ptr), value :: outGradient
      end function cuestXCDerivativeUKSCompute

      integer(c_int) function cuestXCDerivativeUKSComputeWorkspaceQuery( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspaceDescriptor, numOccupiedAlpha, &
         numOccupiedBeta, coefficientMatrixAlpha, coefficientMatrixBeta, outGradient) &
         bind(C, name="cuestXCDerivativeUKSComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int64_t), value :: numOccupiedAlpha
         integer(c_int64_t), value :: numOccupiedBeta
         type(c_ptr), value :: coefficientMatrixAlpha
         type(c_ptr), value :: coefficientMatrixBeta
         type(c_ptr), value :: outGradient
      end function cuestXCDerivativeUKSComputeWorkspaceQuery

      integer(c_int) function cuestXCDerivativeCompute( &
         handle, plan, approximation, parameters, variableBufferSize, temporaryWorkspace, numOccupied, &
         coefficientMatrix, gridXCPotential, outBasisGradient, outGridGradient) &
         bind(C, name="cuestXCDerivativeCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         integer(c_int), value :: approximation
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: coefficientMatrix
         type(c_ptr), value :: gridXCPotential
         type(c_ptr), value :: outBasisGradient
         type(c_ptr), value :: outGridGradient
      end function cuestXCDerivativeCompute

      integer(c_int) function cuestXCDerivativeComputeWorkspaceQuery( &
         handle, plan, approximation, parameters, variableBufferSize, temporaryWorkspaceDescriptor, numOccupied, &
         coefficientMatrix, gridXCPotential, outBasisGradient, outGridGradient) &
         bind(C, name="cuestXCDerivativeComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         integer(c_int), value :: approximation
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: coefficientMatrix
         type(c_ptr), value :: gridXCPotential
         type(c_ptr), value :: outBasisGradient
         type(c_ptr), value :: outGridGradient
      end function cuestXCDerivativeComputeWorkspaceQuery

      integer(c_int) function cuestXCGridDerivativeCompute( &
         handle, plan, parameters, temporaryWorkspace, gridXCEnergy, gridGradient, outGradient) &
         bind(C, name="cuestXCGridDerivativeCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), value :: gridXCEnergy
         type(c_ptr), value :: gridGradient
         type(c_ptr), value :: outGradient
      end function cuestXCGridDerivativeCompute

      integer(c_int) function cuestXCGridDerivativeComputeWorkspaceQuery( &
         handle, plan, parameters, temporaryWorkspaceDescriptor, gridXCEnergy, gridGradient, outGradient) &
         bind(C, name="cuestXCGridDerivativeComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), value :: gridXCEnergy
         type(c_ptr), value :: gridGradient
         type(c_ptr), value :: outGradient
      end function cuestXCGridDerivativeComputeWorkspaceQuery

      integer(c_int) function cuestXCNonsymmetricDerivativeCompute( &
         handle, plan, approximation, parameters, variableBufferSize, temporaryWorkspace, numOccupied, &
         leftCoefficientMatrix, rightCoefficientMatrix, gridXCPotential, outBasisGradient, outGridGradient) &
         bind(C, name="cuestXCNonsymmetricDerivativeCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         integer(c_int), value :: approximation
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: leftCoefficientMatrix
         type(c_ptr), value :: rightCoefficientMatrix
         type(c_ptr), value :: gridXCPotential
         type(c_ptr), value :: outBasisGradient
         type(c_ptr), value :: outGridGradient
      end function cuestXCNonsymmetricDerivativeCompute

      integer(c_int) function cuestXCNonsymmetricDerivativeComputeWorkspaceQuery( &
         handle, plan, approximation, parameters, variableBufferSize, temporaryWorkspaceDescriptor, numOccupied, &
         leftCoefficientMatrix, rightCoefficientMatrix, gridXCPotential, outBasisGradient, outGridGradient) &
         bind(C, name="cuestXCNonsymmetricDerivativeComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         integer(c_int), value :: approximation
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: leftCoefficientMatrix
         type(c_ptr), value :: rightCoefficientMatrix
         type(c_ptr), value :: gridXCPotential
         type(c_ptr), value :: outBasisGradient
         type(c_ptr), value :: outGridGradient
      end function cuestXCNonsymmetricDerivativeComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/exchange_correlation/nonlocal_xc_potential_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestNonlocalXCPotentialRKSCompute( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspace, numOccupied, coefficientMatrix, &
         outXCEnergy, outXCPotentialMatrix) &
         bind(C, name="cuestNonlocalXCPotentialRKSCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: coefficientMatrix
         type(c_ptr), value :: outXCEnergy
         type(c_ptr), value :: outXCPotentialMatrix
      end function cuestNonlocalXCPotentialRKSCompute

      integer(c_int) function cuestNonlocalXCPotentialRKSComputeWorkspaceQuery( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspaceDescriptor, numOccupied, &
         coefficientMatrix, outXCEnergy, outXCPotentialMatrix) &
         bind(C, name="cuestNonlocalXCPotentialRKSComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: coefficientMatrix
         type(c_ptr), value :: outXCEnergy
         type(c_ptr), value :: outXCPotentialMatrix
      end function cuestNonlocalXCPotentialRKSComputeWorkspaceQuery

      integer(c_int) function cuestNonlocalXCPotentialUKSCompute( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspace, numOccupiedAlpha, numOccupiedBeta, &
         coefficientMatrixAlpha, coefficientMatrixBeta, outXCEnergy, outXCPotentialMatrix) &
         bind(C, name="cuestNonlocalXCPotentialUKSCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int64_t), value :: numOccupiedAlpha
         integer(c_int64_t), value :: numOccupiedBeta
         type(c_ptr), value :: coefficientMatrixAlpha
         type(c_ptr), value :: coefficientMatrixBeta
         type(c_ptr), value :: outXCEnergy
         type(c_ptr), value :: outXCPotentialMatrix
      end function cuestNonlocalXCPotentialUKSCompute

      integer(c_int) function cuestNonlocalXCPotentialUKSComputeWorkspaceQuery( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspaceDescriptor, numOccupiedAlpha, &
         numOccupiedBeta, coefficientMatrixAlpha, coefficientMatrixBeta, outXCEnergy, outXCPotentialMatrix) &
         bind(C, name="cuestNonlocalXCPotentialUKSComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int64_t), value :: numOccupiedAlpha
         integer(c_int64_t), value :: numOccupiedBeta
         type(c_ptr), value :: coefficientMatrixAlpha
         type(c_ptr), value :: coefficientMatrixBeta
         type(c_ptr), value :: outXCEnergy
         type(c_ptr), value :: outXCPotentialMatrix
      end function cuestNonlocalXCPotentialUKSComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/exchange_correlation/nonlocal_xc_derivative_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestNonlocalXCDerivativeRKSCompute( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspace, numOccupied, coefficientMatrix, &
         outGradient) &
         bind(C, name="cuestNonlocalXCDerivativeRKSCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: coefficientMatrix
         type(c_ptr), value :: outGradient
      end function cuestNonlocalXCDerivativeRKSCompute

      integer(c_int) function cuestNonlocalXCDerivativeRKSComputeWorkspaceQuery( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspaceDescriptor, numOccupied, &
         coefficientMatrix, outGradient) &
         bind(C, name="cuestNonlocalXCDerivativeRKSComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int64_t), value :: numOccupied
         type(c_ptr), value :: coefficientMatrix
         type(c_ptr), value :: outGradient
      end function cuestNonlocalXCDerivativeRKSComputeWorkspaceQuery

      integer(c_int) function cuestNonlocalXCDerivativeUKSCompute( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspace, numOccupiedAlpha, numOccupiedBeta, &
         coefficientMatrixAlpha, coefficientMatrixBeta, outGradient) &
         bind(C, name="cuestNonlocalXCDerivativeUKSCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int64_t), value :: numOccupiedAlpha
         integer(c_int64_t), value :: numOccupiedBeta
         type(c_ptr), value :: coefficientMatrixAlpha
         type(c_ptr), value :: coefficientMatrixBeta
         type(c_ptr), value :: outGradient
      end function cuestNonlocalXCDerivativeUKSCompute

      integer(c_int) function cuestNonlocalXCDerivativeUKSComputeWorkspaceQuery( &
         handle, plan, parameters, variableBufferSize, temporaryWorkspaceDescriptor, numOccupiedAlpha, &
         numOccupiedBeta, coefficientMatrixAlpha, coefficientMatrixBeta, outGradient) &
         bind(C, name="cuestNonlocalXCDerivativeUKSComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(in) :: variableBufferSize
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int64_t), value :: numOccupiedAlpha
         integer(c_int64_t), value :: numOccupiedBeta
         type(c_ptr), value :: coefficientMatrixAlpha
         type(c_ptr), value :: coefficientMatrixBeta
         type(c_ptr), value :: outGradient
      end function cuestNonlocalXCDerivativeUKSComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/integral_plans/pcm_integral_plan_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestPCMIntPlanCreate( &
         handle, intPlan, parameters, persistentWorkspace, temporaryWorkspace, numAngularPointsPerAtom, epsilon, &
         zetas, atomicRadii, effectiveNuclearCharges, outPlan) &
         bind(C, name="cuestPCMIntPlanCreate")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: intPlan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: persistentWorkspace
         type(cuestWorkspace_t) :: temporaryWorkspace
         integer(c_int64_t), dimension(*), intent(in) :: numAngularPointsPerAtom
         real(c_double), value :: epsilon
         type(c_ptr), value :: zetas
         type(c_ptr), value :: atomicRadii
         type(c_ptr), value :: effectiveNuclearCharges
         type(c_ptr), intent(out) :: outPlan
      end function cuestPCMIntPlanCreate

      integer(c_int) function cuestPCMIntPlanCreateWorkspaceQuery( &
         handle, intPlan, parameters, persistentWorkspaceDescriptor, temporaryWorkspaceDescriptor, &
         numAngularPointsPerAtom, epsilon, zetas, atomicRadii, effectiveNuclearCharges, outPlan) &
         bind(C, name="cuestPCMIntPlanCreateWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: intPlan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: persistentWorkspaceDescriptor
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         integer(c_int64_t), dimension(*), intent(in) :: numAngularPointsPerAtom
         real(c_double), value :: epsilon
         type(c_ptr), value :: zetas
         type(c_ptr), value :: atomicRadii
         type(c_ptr), value :: effectiveNuclearCharges
         type(c_ptr), intent(out) :: outPlan
      end function cuestPCMIntPlanCreateWorkspaceQuery

      integer(c_int) function cuestPCMIntPlanDestroy(plan) &
         bind(C, name="cuestPCMIntPlanDestroy")
         import
         type(c_ptr), value :: plan
      end function cuestPCMIntPlanDestroy

      ! ---------------------------------------------------------------
      !  cuest/pcm/pcm_potential_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestPCMPotentialCompute( &
         handle, plan, parameters, temporaryWorkspace, densityMatrix, inQ, outQ, outPCMResults, &
         outPCMPotentialMatrix) &
         bind(C, name="cuestPCMPotentialCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: inQ
         type(c_ptr), value :: outQ
         type(c_ptr), value :: outPCMResults
         type(c_ptr), value :: outPCMPotentialMatrix
      end function cuestPCMPotentialCompute

      integer(c_int) function cuestPCMPotentialComputeWorkspaceQuery( &
         handle, plan, parameters, temporaryWorkspaceDescriptor, densityMatrix, inQ, outQ, outPCMResults, &
         outPCMPotentialMatrix) &
         bind(C, name="cuestPCMPotentialComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: inQ
         type(c_ptr), value :: outQ
         type(c_ptr), value :: outPCMResults
         type(c_ptr), value :: outPCMPotentialMatrix
      end function cuestPCMPotentialComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/pcm/pcm_derivative_compute_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestPCMDerivativeCompute( &
         handle, plan, parameters, temporaryWorkspace, densityMatrix, inQ, outQ, outPCMResults, outPCMGradient) &
         bind(C, name="cuestPCMDerivativeCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: inQ
         type(c_ptr), value :: outQ
         type(c_ptr), value :: outPCMResults
         type(c_ptr), value :: outPCMGradient
      end function cuestPCMDerivativeCompute

      integer(c_int) function cuestPCMDerivativeComputeWorkspaceQuery( &
         handle, plan, parameters, temporaryWorkspaceDescriptor, densityMatrix, inQ, outQ, outPCMResults, &
         outPCMGradient) &
         bind(C, name="cuestPCMDerivativeComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: inQ
         type(c_ptr), value :: outQ
         type(c_ptr), value :: outPCMResults
         type(c_ptr), value :: outPCMGradient
      end function cuestPCMDerivativeComputeWorkspaceQuery

      integer(c_int) function cuestPCMRadiiDerivativeCompute( &
         handle, plan, parameters, temporaryWorkspace, densityMatrix, inQ, outQ, outPCMResults, outPCMRadiiGradient) &
         bind(C, name="cuestPCMRadiiDerivativeCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: inQ
         type(c_ptr), value :: outQ
         type(c_ptr), value :: outPCMResults
         type(c_ptr), value :: outPCMRadiiGradient
      end function cuestPCMRadiiDerivativeCompute

      integer(c_int) function cuestPCMRadiiDerivativeComputeWorkspaceQuery( &
         handle, plan, parameters, temporaryWorkspaceDescriptor, densityMatrix, inQ, outQ, outPCMResults, &
         outPCMRadiiGradient) &
         bind(C, name="cuestPCMRadiiDerivativeComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), value :: densityMatrix
         type(c_ptr), value :: inQ
         type(c_ptr), value :: outQ
         type(c_ptr), value :: outPCMResults
         type(c_ptr), value :: outPCMRadiiGradient
      end function cuestPCMRadiiDerivativeComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/pcm/pcm_integration_grid_compute_api.h
      ! ---------------------------------------------------------------
     integer(c_int) function cuestPCMIntegrationGridCompute(handle, plan, parameters, temporaryWorkspace, outGridPoints) &
         bind(C, name="cuestPCMIntegrationGridCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), value :: outGridPoints
      end function cuestPCMIntegrationGridCompute

      integer(c_int) function cuestPCMIntegrationGridComputeWorkspaceQuery( &
         handle, plan, parameters, temporaryWorkspaceDescriptor, outGridPoints) &
         bind(C, name="cuestPCMIntegrationGridComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), value :: outGridPoints
      end function cuestPCMIntegrationGridComputeWorkspaceQuery

      integer(c_int) function cuestPCMIntegrationWeightCompute( &
         handle, plan, weightType, parameters, temporaryWorkspace, outGridWeights) &
         bind(C, name="cuestPCMIntegrationWeightCompute")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         integer(c_int), value :: weightType
         type(c_ptr), value :: parameters
         type(cuestWorkspace_t) :: temporaryWorkspace
         type(c_ptr), value :: outGridWeights
      end function cuestPCMIntegrationWeightCompute

      integer(c_int) function cuestPCMIntegrationWeightComputeWorkspaceQuery( &
         handle, plan, weightType, parameters, temporaryWorkspaceDescriptor, outGridWeights) &
         bind(C, name="cuestPCMIntegrationWeightComputeWorkspaceQuery")
         import
         type(c_ptr), value :: handle
         type(c_ptr), value :: plan
         integer(c_int), value :: weightType
         type(c_ptr), value :: parameters
         type(cuestWorkspaceDescriptor_t), intent(out) :: temporaryWorkspaceDescriptor
         type(c_ptr), value :: outGridWeights
      end function cuestPCMIntegrationWeightComputeWorkspaceQuery

      ! ---------------------------------------------------------------
      !  cuest/query/query_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestQuery(handle, type, object, attribute, attributeValue, attributeValueSize) &
         bind(C, name="cuestQuery")
         import
         type(c_ptr), value :: handle
         integer(c_int), value :: type
         type(c_ptr), value :: object
         integer(c_int), value :: attribute
         type(c_ptr), value :: attributeValue
         integer(c_size_t), value :: attributeValueSize
      end function cuestQuery

      ! ---------------------------------------------------------------
      !  cuest/results/results_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestResultsCreate(resultsType, outResults) &
         bind(C, name="cuestResultsCreate")
         import
         integer(c_int), value :: resultsType
         type(c_ptr), intent(out) :: outResults
      end function cuestResultsCreate

      integer(c_int) function cuestResultsDestroy(resultsType, results) &
         bind(C, name="cuestResultsDestroy")
         import
         integer(c_int), value :: resultsType
         type(c_ptr), value :: results
      end function cuestResultsDestroy

      integer(c_int) function cuestResultsQuery(resultsType, results, attribute, resultValue, resultValueSize) &
         bind(C, name="cuestResultsQuery")
         import
         integer(c_int), value :: resultsType
         type(c_ptr), value :: results
         integer(c_int), value :: attribute
         type(c_ptr), value :: resultValue
         integer(c_size_t), value :: resultValueSize
      end function cuestResultsQuery

      ! ---------------------------------------------------------------
      !  cuest/parameters/parameter_api.h
      ! ---------------------------------------------------------------
      integer(c_int) function cuestParametersCreate(parameterType, outParameters) &
         bind(C, name="cuestParametersCreate")
         import
         integer(c_int), value :: parameterType
         type(c_ptr), intent(out) :: outParameters
      end function cuestParametersCreate

      integer(c_int) function cuestParametersDestroy(parameterType, parameters) &
         bind(C, name="cuestParametersDestroy")
         import
         integer(c_int), value :: parameterType
         type(c_ptr), value :: parameters
      end function cuestParametersDestroy

  integer(c_int) function cuestParametersQuery(parameterType, parameters, attribute, attributeValue, attributeValueSize) &
         bind(C, name="cuestParametersQuery")
         import
         integer(c_int), value :: parameterType
         type(c_ptr), value :: parameters
         integer(c_int), value :: attribute
         type(c_ptr), value :: attributeValue
         integer(c_size_t), value :: attributeValueSize
      end function cuestParametersQuery

      integer(c_int) function cuestParametersConfigure( &
         parameterType, parameters, attribute, attributeValue, attributeValueSize) &
         bind(C, name="cuestParametersConfigure")
         import
         integer(c_int), value :: parameterType
         type(c_ptr), value :: parameters
         integer(c_int), value :: attribute
         type(c_ptr), value :: attributeValue
         integer(c_size_t), value :: attributeValueSize
      end function cuestParametersConfigure

   end interface

end module cuest
