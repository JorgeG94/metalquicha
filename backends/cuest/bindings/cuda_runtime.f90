! ==========================================================================
!  cuda_runtime.f90 -- Fortran 2008 iso_c_binding interface to the CUDA Runtime API
!
!  GENERATED FILE -- do not edit by hand.
!  Regenerate with:  python3 generate_cuda_fortran.py
!  Generated from CUDA 12.9 headers at:
!      /apps/cuda/12.9.0/include
!
!  Standard Fortran 2008 only -- no compiler extensions. Builds with
!  gfortran, ifx, flang and nvfortran, and is independent of cudafor.
!
!  Conventions
!  -----------
!  * Every entry point returns its C status code as an INTEGER(c_int)
!    function result (cudaError_t / CUresult).
!  * Opaque handles (streams, events, graphs, arrays, ...) are pointers
!    in C and map to TYPE(c_ptr): by VALUE when passed in, INTENT(OUT)
!    when returned (C  T*  where T is itself a pointer typedef).
!  * Device and host buffers (void*, T*) are TYPE(c_ptr), VALUE. Pass a
!    device address, or C_LOC(host_array) for host data.
!  * Enumerators are PUBLIC INTEGER(c_int) PARAMETERs with their C names.
!  * Structs are BIND(C) derived types whose layout has been verified
!    against sizeof/offsetof from a compiled C probe. C unions have no
!    Fortran equivalent and appear as INTEGER(c_int8_t) byte arrays of
!    the correct size.
!  * Fortran is case-insensitive; the C names are preserved verbatim and
!    checked for case-insensitive collisions at generation time.
! ==========================================================================
module cuda_runtime
   use, intrinsic :: iso_c_binding
   implicit none
   public

   ! ======================================================================
   !  Enumerations
   ! ======================================================================
   ! ---- cudaRoundMode
   integer(c_int), parameter :: cudaRoundNearest = 0
   integer(c_int), parameter :: cudaRoundZero = 1
   integer(c_int), parameter :: cudaRoundPosInf = 2
   integer(c_int), parameter :: cudaRoundMinInf = 3

   ! ---- cudaError
   integer(c_int), parameter :: cudaSuccess = 0
   integer(c_int), parameter :: cudaErrorInvalidValue = 1
   integer(c_int), parameter :: cudaErrorMemoryAllocation = 2
   integer(c_int), parameter :: cudaErrorInitializationError = 3
   integer(c_int), parameter :: cudaErrorCudartUnloading = 4
   integer(c_int), parameter :: cudaErrorProfilerDisabled = 5
   integer(c_int), parameter :: cudaErrorProfilerNotInitialized = 6
   integer(c_int), parameter :: cudaErrorProfilerAlreadyStarted = 7
   integer(c_int), parameter :: cudaErrorProfilerAlreadyStopped = 8
   integer(c_int), parameter :: cudaErrorInvalidConfiguration = 9
   integer(c_int), parameter :: cudaErrorInvalidPitchValue = 12
   integer(c_int), parameter :: cudaErrorInvalidSymbol = 13
   integer(c_int), parameter :: cudaErrorInvalidHostPointer = 16
   integer(c_int), parameter :: cudaErrorInvalidDevicePointer = 17
   integer(c_int), parameter :: cudaErrorInvalidTexture = 18
   integer(c_int), parameter :: cudaErrorInvalidTextureBinding = 19
   integer(c_int), parameter :: cudaErrorInvalidChannelDescriptor = 20
   integer(c_int), parameter :: cudaErrorInvalidMemcpyDirection = 21
   integer(c_int), parameter :: cudaErrorAddressOfConstant = 22
   integer(c_int), parameter :: cudaErrorTextureFetchFailed = 23
   integer(c_int), parameter :: cudaErrorTextureNotBound = 24
   integer(c_int), parameter :: cudaErrorSynchronizationError = 25
   integer(c_int), parameter :: cudaErrorInvalidFilterSetting = 26
   integer(c_int), parameter :: cudaErrorInvalidNormSetting = 27
   integer(c_int), parameter :: cudaErrorMixedDeviceExecution = 28
   integer(c_int), parameter :: cudaErrorNotYetImplemented = 31
   integer(c_int), parameter :: cudaErrorMemoryValueTooLarge = 32
   integer(c_int), parameter :: cudaErrorStubLibrary = 34
   integer(c_int), parameter :: cudaErrorInsufficientDriver = 35
   integer(c_int), parameter :: cudaErrorCallRequiresNewerDriver = 36
   integer(c_int), parameter :: cudaErrorInvalidSurface = 37
   integer(c_int), parameter :: cudaErrorDuplicateVariableName = 43
   integer(c_int), parameter :: cudaErrorDuplicateTextureName = 44
   integer(c_int), parameter :: cudaErrorDuplicateSurfaceName = 45
   integer(c_int), parameter :: cudaErrorDevicesUnavailable = 46
   integer(c_int), parameter :: cudaErrorIncompatibleDriverContext = 49
   integer(c_int), parameter :: cudaErrorMissingConfiguration = 52
   integer(c_int), parameter :: cudaErrorPriorLaunchFailure = 53
   integer(c_int), parameter :: cudaErrorLaunchMaxDepthExceeded = 65
   integer(c_int), parameter :: cudaErrorLaunchFileScopedTex = 66
   integer(c_int), parameter :: cudaErrorLaunchFileScopedSurf = 67
   integer(c_int), parameter :: cudaErrorSyncDepthExceeded = 68
   integer(c_int), parameter :: cudaErrorLaunchPendingCountExceeded = 69
   integer(c_int), parameter :: cudaErrorInvalidDeviceFunction = 98
   integer(c_int), parameter :: cudaErrorNoDevice = 100
   integer(c_int), parameter :: cudaErrorInvalidDevice = 101
   integer(c_int), parameter :: cudaErrorDeviceNotLicensed = 102
   integer(c_int), parameter :: cudaErrorSoftwareValidityNotEstablished = 103
   integer(c_int), parameter :: cudaErrorStartupFailure = 127
   integer(c_int), parameter :: cudaErrorInvalidKernelImage = 200
   integer(c_int), parameter :: cudaErrorDeviceUninitialized = 201
   integer(c_int), parameter :: cudaErrorMapBufferObjectFailed = 205
   integer(c_int), parameter :: cudaErrorUnmapBufferObjectFailed = 206
   integer(c_int), parameter :: cudaErrorArrayIsMapped = 207
   integer(c_int), parameter :: cudaErrorAlreadyMapped = 208
   integer(c_int), parameter :: cudaErrorNoKernelImageForDevice = 209
   integer(c_int), parameter :: cudaErrorAlreadyAcquired = 210
   integer(c_int), parameter :: cudaErrorNotMapped = 211
   integer(c_int), parameter :: cudaErrorNotMappedAsArray = 212
   integer(c_int), parameter :: cudaErrorNotMappedAsPointer = 213
   integer(c_int), parameter :: cudaErrorECCUncorrectable = 214
   integer(c_int), parameter :: cudaErrorUnsupportedLimit = 215
   integer(c_int), parameter :: cudaErrorDeviceAlreadyInUse = 216
   integer(c_int), parameter :: cudaErrorPeerAccessUnsupported = 217
   integer(c_int), parameter :: cudaErrorInvalidPtx = 218
   integer(c_int), parameter :: cudaErrorInvalidGraphicsContext = 219
   integer(c_int), parameter :: cudaErrorNvlinkUncorrectable = 220
   integer(c_int), parameter :: cudaErrorJitCompilerNotFound = 221
   integer(c_int), parameter :: cudaErrorUnsupportedPtxVersion = 222
   integer(c_int), parameter :: cudaErrorJitCompilationDisabled = 223
   integer(c_int), parameter :: cudaErrorUnsupportedExecAffinity = 224
   integer(c_int), parameter :: cudaErrorUnsupportedDevSideSync = 225
   integer(c_int), parameter :: cudaErrorContained = 226
   integer(c_int), parameter :: cudaErrorInvalidSource = 300
   integer(c_int), parameter :: cudaErrorFileNotFound = 301
   integer(c_int), parameter :: cudaErrorSharedObjectSymbolNotFound = 302
   integer(c_int), parameter :: cudaErrorSharedObjectInitFailed = 303
   integer(c_int), parameter :: cudaErrorOperatingSystem = 304
   integer(c_int), parameter :: cudaErrorInvalidResourceHandle = 400
   integer(c_int), parameter :: cudaErrorIllegalState = 401
   integer(c_int), parameter :: cudaErrorLossyQuery = 402
   integer(c_int), parameter :: cudaErrorSymbolNotFound = 500
   integer(c_int), parameter :: cudaErrorNotReady = 600
   integer(c_int), parameter :: cudaErrorIllegalAddress = 700
   integer(c_int), parameter :: cudaErrorLaunchOutOfResources = 701
   integer(c_int), parameter :: cudaErrorLaunchTimeout = 702
   integer(c_int), parameter :: cudaErrorLaunchIncompatibleTexturing = 703
   integer(c_int), parameter :: cudaErrorPeerAccessAlreadyEnabled = 704
   integer(c_int), parameter :: cudaErrorPeerAccessNotEnabled = 705
   integer(c_int), parameter :: cudaErrorSetOnActiveProcess = 708
   integer(c_int), parameter :: cudaErrorContextIsDestroyed = 709
   integer(c_int), parameter :: cudaErrorAssert = 710
   integer(c_int), parameter :: cudaErrorTooManyPeers = 711
   integer(c_int), parameter :: cudaErrorHostMemoryAlreadyRegistered = 712
   integer(c_int), parameter :: cudaErrorHostMemoryNotRegistered = 713
   integer(c_int), parameter :: cudaErrorHardwareStackError = 714
   integer(c_int), parameter :: cudaErrorIllegalInstruction = 715
   integer(c_int), parameter :: cudaErrorMisalignedAddress = 716
   integer(c_int), parameter :: cudaErrorInvalidAddressSpace = 717
   integer(c_int), parameter :: cudaErrorInvalidPc = 718
   integer(c_int), parameter :: cudaErrorLaunchFailure = 719
   integer(c_int), parameter :: cudaErrorCooperativeLaunchTooLarge = 720
   integer(c_int), parameter :: cudaErrorTensorMemoryLeak = 721
   integer(c_int), parameter :: cudaErrorNotPermitted = 800
   integer(c_int), parameter :: cudaErrorNotSupported = 801
   integer(c_int), parameter :: cudaErrorSystemNotReady = 802
   integer(c_int), parameter :: cudaErrorSystemDriverMismatch = 803
   integer(c_int), parameter :: cudaErrorCompatNotSupportedOnDevice = 804
   integer(c_int), parameter :: cudaErrorMpsConnectionFailed = 805
   integer(c_int), parameter :: cudaErrorMpsRpcFailure = 806
   integer(c_int), parameter :: cudaErrorMpsServerNotReady = 807
   integer(c_int), parameter :: cudaErrorMpsMaxClientsReached = 808
   integer(c_int), parameter :: cudaErrorMpsMaxConnectionsReached = 809
   integer(c_int), parameter :: cudaErrorMpsClientTerminated = 810
   integer(c_int), parameter :: cudaErrorCdpNotSupported = 811
   integer(c_int), parameter :: cudaErrorCdpVersionMismatch = 812
   integer(c_int), parameter :: cudaErrorStreamCaptureUnsupported = 900
   integer(c_int), parameter :: cudaErrorStreamCaptureInvalidated = 901
   integer(c_int), parameter :: cudaErrorStreamCaptureMerge = 902
   integer(c_int), parameter :: cudaErrorStreamCaptureUnmatched = 903
   integer(c_int), parameter :: cudaErrorStreamCaptureUnjoined = 904
   integer(c_int), parameter :: cudaErrorStreamCaptureIsolation = 905
   integer(c_int), parameter :: cudaErrorStreamCaptureImplicit = 906
   integer(c_int), parameter :: cudaErrorCapturedEvent = 907
   integer(c_int), parameter :: cudaErrorStreamCaptureWrongThread = 908
   integer(c_int), parameter :: cudaErrorTimeout = 909
   integer(c_int), parameter :: cudaErrorGraphExecUpdateFailure = 910
   integer(c_int), parameter :: cudaErrorExternalDevice = 911
   integer(c_int), parameter :: cudaErrorInvalidClusterSize = 912
   integer(c_int), parameter :: cudaErrorFunctionNotLoaded = 913
   integer(c_int), parameter :: cudaErrorInvalidResourceType = 914
   integer(c_int), parameter :: cudaErrorInvalidResourceConfiguration = 915
   integer(c_int), parameter :: cudaErrorUnknown = 999
   integer(c_int), parameter :: cudaErrorApiFailureBase = 10000

   ! ---- cudaChannelFormatKind
   integer(c_int), parameter :: cudaChannelFormatKindSigned = 0
   integer(c_int), parameter :: cudaChannelFormatKindUnsigned = 1
   integer(c_int), parameter :: cudaChannelFormatKindFloat = 2
   integer(c_int), parameter :: cudaChannelFormatKindNone = 3
   integer(c_int), parameter :: cudaChannelFormatKindNV12 = 4
   integer(c_int), parameter :: cudaChannelFormatKindUnsignedNormalized8X1 = 5
   integer(c_int), parameter :: cudaChannelFormatKindUnsignedNormalized8X2 = 6
   integer(c_int), parameter :: cudaChannelFormatKindUnsignedNormalized8X4 = 7
   integer(c_int), parameter :: cudaChannelFormatKindUnsignedNormalized16X1 = 8
   integer(c_int), parameter :: cudaChannelFormatKindUnsignedNormalized16X2 = 9
   integer(c_int), parameter :: cudaChannelFormatKindUnsignedNormalized16X4 = 10
   integer(c_int), parameter :: cudaChannelFormatKindSignedNormalized8X1 = 11
   integer(c_int), parameter :: cudaChannelFormatKindSignedNormalized8X2 = 12
   integer(c_int), parameter :: cudaChannelFormatKindSignedNormalized8X4 = 13
   integer(c_int), parameter :: cudaChannelFormatKindSignedNormalized16X1 = 14
   integer(c_int), parameter :: cudaChannelFormatKindSignedNormalized16X2 = 15
   integer(c_int), parameter :: cudaChannelFormatKindSignedNormalized16X4 = 16
   integer(c_int), parameter :: cudaChannelFormatKindUnsignedBlockCompressed1 = 17
   integer(c_int), parameter :: cudaChannelFormatKindUnsignedBlockCompressed1SRGB = 18
   integer(c_int), parameter :: cudaChannelFormatKindUnsignedBlockCompressed2 = 19
   integer(c_int), parameter :: cudaChannelFormatKindUnsignedBlockCompressed2SRGB = 20
   integer(c_int), parameter :: cudaChannelFormatKindUnsignedBlockCompressed3 = 21
   integer(c_int), parameter :: cudaChannelFormatKindUnsignedBlockCompressed3SRGB = 22
   integer(c_int), parameter :: cudaChannelFormatKindUnsignedBlockCompressed4 = 23
   integer(c_int), parameter :: cudaChannelFormatKindSignedBlockCompressed4 = 24
   integer(c_int), parameter :: cudaChannelFormatKindUnsignedBlockCompressed5 = 25
   integer(c_int), parameter :: cudaChannelFormatKindSignedBlockCompressed5 = 26
   integer(c_int), parameter :: cudaChannelFormatKindUnsignedBlockCompressed6H = 27
   integer(c_int), parameter :: cudaChannelFormatKindSignedBlockCompressed6H = 28
   integer(c_int), parameter :: cudaChannelFormatKindUnsignedBlockCompressed7 = 29
   integer(c_int), parameter :: cudaChannelFormatKindUnsignedBlockCompressed7SRGB = 30
   integer(c_int), parameter :: cudaChannelFormatKindUnsignedNormalized1010102 = 31

   ! ---- cudaMemoryType
   integer(c_int), parameter :: cudaMemoryTypeUnregistered = 0
   integer(c_int), parameter :: cudaMemoryTypeHost = 1
   integer(c_int), parameter :: cudaMemoryTypeDevice = 2
   integer(c_int), parameter :: cudaMemoryTypeManaged = 3

   ! ---- cudaMemcpyKind
   integer(c_int), parameter :: cudaMemcpyHostToHost = 0
   integer(c_int), parameter :: cudaMemcpyHostToDevice = 1
   integer(c_int), parameter :: cudaMemcpyDeviceToHost = 2
   integer(c_int), parameter :: cudaMemcpyDeviceToDevice = 3
   integer(c_int), parameter :: cudaMemcpyDefault = 4

   ! ---- cudaAccessProperty
   integer(c_int), parameter :: cudaAccessPropertyNormal = 0
   integer(c_int), parameter :: cudaAccessPropertyStreaming = 1
   integer(c_int), parameter :: cudaAccessPropertyPersisting = 2

   ! ---- cudaStreamCaptureStatus
   integer(c_int), parameter :: cudaStreamCaptureStatusNone = 0
   integer(c_int), parameter :: cudaStreamCaptureStatusActive = 1
   integer(c_int), parameter :: cudaStreamCaptureStatusInvalidated = 2

   ! ---- cudaStreamCaptureMode
   integer(c_int), parameter :: cudaStreamCaptureModeGlobal = 0
   integer(c_int), parameter :: cudaStreamCaptureModeThreadLocal = 1
   integer(c_int), parameter :: cudaStreamCaptureModeRelaxed = 2

   ! ---- cudaSynchronizationPolicy
   integer(c_int), parameter :: cudaSyncPolicyAuto = 1
   integer(c_int), parameter :: cudaSyncPolicySpin = 2
   integer(c_int), parameter :: cudaSyncPolicyYield = 3
   integer(c_int), parameter :: cudaSyncPolicyBlockingSync = 4

   ! ---- cudaClusterSchedulingPolicy
   integer(c_int), parameter :: cudaClusterSchedulingPolicyDefault = 0
   integer(c_int), parameter :: cudaClusterSchedulingPolicySpread = 1
   integer(c_int), parameter :: cudaClusterSchedulingPolicyLoadBalancing = 2

   ! ---- cudaStreamUpdateCaptureDependenciesFlags
   integer(c_int), parameter :: cudaStreamAddCaptureDependencies = 0
   integer(c_int), parameter :: cudaStreamSetCaptureDependencies = 1

   ! ---- cudaUserObjectFlags
   integer(c_int), parameter :: cudaUserObjectNoDestructorSync = 1

   ! ---- cudaUserObjectRetainFlags
   integer(c_int), parameter :: cudaGraphUserObjectMove = 1

   ! ---- cudaGraphicsRegisterFlags
   integer(c_int), parameter :: cudaGraphicsRegisterFlagsNone = 0
   integer(c_int), parameter :: cudaGraphicsRegisterFlagsReadOnly = 1
   integer(c_int), parameter :: cudaGraphicsRegisterFlagsWriteDiscard = 2
   integer(c_int), parameter :: cudaGraphicsRegisterFlagsSurfaceLoadStore = 4
   integer(c_int), parameter :: cudaGraphicsRegisterFlagsTextureGather = 8

   ! ---- cudaGraphicsMapFlags
   integer(c_int), parameter :: cudaGraphicsMapFlagsNone = 0
   integer(c_int), parameter :: cudaGraphicsMapFlagsReadOnly = 1
   integer(c_int), parameter :: cudaGraphicsMapFlagsWriteDiscard = 2

   ! ---- cudaGraphicsCubeFace
   integer(c_int), parameter :: cudaGraphicsCubeFacePositiveX = 0
   integer(c_int), parameter :: cudaGraphicsCubeFaceNegativeX = 1
   integer(c_int), parameter :: cudaGraphicsCubeFacePositiveY = 2
   integer(c_int), parameter :: cudaGraphicsCubeFaceNegativeY = 3
   integer(c_int), parameter :: cudaGraphicsCubeFacePositiveZ = 4
   integer(c_int), parameter :: cudaGraphicsCubeFaceNegativeZ = 5

   ! ---- cudaResourceType
   integer(c_int), parameter :: cudaResourceTypeArray = 0
   integer(c_int), parameter :: cudaResourceTypeMipmappedArray = 1
   integer(c_int), parameter :: cudaResourceTypeLinear = 2
   integer(c_int), parameter :: cudaResourceTypePitch2D = 3

   ! ---- cudaResourceViewFormat
   integer(c_int), parameter :: cudaResViewFormatNone = 0
   integer(c_int), parameter :: cudaResViewFormatUnsignedChar1 = 1
   integer(c_int), parameter :: cudaResViewFormatUnsignedChar2 = 2
   integer(c_int), parameter :: cudaResViewFormatUnsignedChar4 = 3
   integer(c_int), parameter :: cudaResViewFormatSignedChar1 = 4
   integer(c_int), parameter :: cudaResViewFormatSignedChar2 = 5
   integer(c_int), parameter :: cudaResViewFormatSignedChar4 = 6
   integer(c_int), parameter :: cudaResViewFormatUnsignedShort1 = 7
   integer(c_int), parameter :: cudaResViewFormatUnsignedShort2 = 8
   integer(c_int), parameter :: cudaResViewFormatUnsignedShort4 = 9
   integer(c_int), parameter :: cudaResViewFormatSignedShort1 = 10
   integer(c_int), parameter :: cudaResViewFormatSignedShort2 = 11
   integer(c_int), parameter :: cudaResViewFormatSignedShort4 = 12
   integer(c_int), parameter :: cudaResViewFormatUnsignedInt1 = 13
   integer(c_int), parameter :: cudaResViewFormatUnsignedInt2 = 14
   integer(c_int), parameter :: cudaResViewFormatUnsignedInt4 = 15
   integer(c_int), parameter :: cudaResViewFormatSignedInt1 = 16
   integer(c_int), parameter :: cudaResViewFormatSignedInt2 = 17
   integer(c_int), parameter :: cudaResViewFormatSignedInt4 = 18
   integer(c_int), parameter :: cudaResViewFormatHalf1 = 19
   integer(c_int), parameter :: cudaResViewFormatHalf2 = 20
   integer(c_int), parameter :: cudaResViewFormatHalf4 = 21
   integer(c_int), parameter :: cudaResViewFormatFloat1 = 22
   integer(c_int), parameter :: cudaResViewFormatFloat2 = 23
   integer(c_int), parameter :: cudaResViewFormatFloat4 = 24
   integer(c_int), parameter :: cudaResViewFormatUnsignedBlockCompressed1 = 25
   integer(c_int), parameter :: cudaResViewFormatUnsignedBlockCompressed2 = 26
   integer(c_int), parameter :: cudaResViewFormatUnsignedBlockCompressed3 = 27
   integer(c_int), parameter :: cudaResViewFormatUnsignedBlockCompressed4 = 28
   integer(c_int), parameter :: cudaResViewFormatSignedBlockCompressed4 = 29
   integer(c_int), parameter :: cudaResViewFormatUnsignedBlockCompressed5 = 30
   integer(c_int), parameter :: cudaResViewFormatSignedBlockCompressed5 = 31
   integer(c_int), parameter :: cudaResViewFormatUnsignedBlockCompressed6H = 32
   integer(c_int), parameter :: cudaResViewFormatSignedBlockCompressed6H = 33
   integer(c_int), parameter :: cudaResViewFormatUnsignedBlockCompressed7 = 34

   ! ---- cudaFuncAttribute
   integer(c_int), parameter :: cudaFuncAttributeMaxDynamicSharedMemorySize = 8
   integer(c_int), parameter :: cudaFuncAttributePreferredSharedMemoryCarveout = 9
   integer(c_int), parameter :: cudaFuncAttributeClusterDimMustBeSet = 10
   integer(c_int), parameter :: cudaFuncAttributeRequiredClusterWidth = 11
   integer(c_int), parameter :: cudaFuncAttributeRequiredClusterHeight = 12
   integer(c_int), parameter :: cudaFuncAttributeRequiredClusterDepth = 13
   integer(c_int), parameter :: cudaFuncAttributeNonPortableClusterSizeAllowed = 14
   integer(c_int), parameter :: cudaFuncAttributeClusterSchedulingPolicyPreference = 15
   integer(c_int), parameter :: cudaFuncAttributeMax = 16

   ! ---- cudaFuncCache
   integer(c_int), parameter :: cudaFuncCachePreferNone = 0
   integer(c_int), parameter :: cudaFuncCachePreferShared = 1
   integer(c_int), parameter :: cudaFuncCachePreferL1 = 2
   integer(c_int), parameter :: cudaFuncCachePreferEqual = 3

   ! ---- cudaSharedMemConfig
   integer(c_int), parameter :: cudaSharedMemBankSizeDefault = 0
   integer(c_int), parameter :: cudaSharedMemBankSizeFourByte = 1
   integer(c_int), parameter :: cudaSharedMemBankSizeEightByte = 2

   ! ---- cudaSharedCarveout
   integer(c_int), parameter :: cudaSharedmemCarveoutDefault = -1
   integer(c_int), parameter :: cudaSharedmemCarveoutMaxShared = 100
   integer(c_int), parameter :: cudaSharedmemCarveoutMaxL1 = 0

   ! ---- cudaComputeMode
   integer(c_int), parameter :: cudaComputeModeDefault = 0
   integer(c_int), parameter :: cudaComputeModeExclusive = 1
   integer(c_int), parameter :: cudaComputeModeProhibited = 2
   integer(c_int), parameter :: cudaComputeModeExclusiveProcess = 3

   ! ---- cudaLimit
   integer(c_int), parameter :: cudaLimitStackSize = 0
   integer(c_int), parameter :: cudaLimitPrintfFifoSize = 1
   integer(c_int), parameter :: cudaLimitMallocHeapSize = 2
   integer(c_int), parameter :: cudaLimitDevRuntimeSyncDepth = 3
   integer(c_int), parameter :: cudaLimitDevRuntimePendingLaunchCount = 4
   integer(c_int), parameter :: cudaLimitMaxL2FetchGranularity = 5
   integer(c_int), parameter :: cudaLimitPersistingL2CacheSize = 6

   ! ---- cudaMemoryAdvise
   integer(c_int), parameter :: cudaMemAdviseSetReadMostly = 1
   integer(c_int), parameter :: cudaMemAdviseUnsetReadMostly = 2
   integer(c_int), parameter :: cudaMemAdviseSetPreferredLocation = 3
   integer(c_int), parameter :: cudaMemAdviseUnsetPreferredLocation = 4
   integer(c_int), parameter :: cudaMemAdviseSetAccessedBy = 5
   integer(c_int), parameter :: cudaMemAdviseUnsetAccessedBy = 6

   ! ---- cudaMemRangeAttribute
   integer(c_int), parameter :: cudaMemRangeAttributeReadMostly = 1
   integer(c_int), parameter :: cudaMemRangeAttributePreferredLocation = 2
   integer(c_int), parameter :: cudaMemRangeAttributeAccessedBy = 3
   integer(c_int), parameter :: cudaMemRangeAttributeLastPrefetchLocation = 4
   integer(c_int), parameter :: cudaMemRangeAttributePreferredLocationType = 5
   integer(c_int), parameter :: cudaMemRangeAttributePreferredLocationId = 6
   integer(c_int), parameter :: cudaMemRangeAttributeLastPrefetchLocationType = 7
   integer(c_int), parameter :: cudaMemRangeAttributeLastPrefetchLocationId = 8

   ! ---- cudaFlushGPUDirectRDMAWritesOptions
   integer(c_int), parameter :: cudaFlushGPUDirectRDMAWritesOptionHost = 1
   integer(c_int), parameter :: cudaFlushGPUDirectRDMAWritesOptionMemOps = 2

   ! ---- cudaGPUDirectRDMAWritesOrdering
   integer(c_int), parameter :: cudaGPUDirectRDMAWritesOrderingNone = 0
   integer(c_int), parameter :: cudaGPUDirectRDMAWritesOrderingOwner = 100
   integer(c_int), parameter :: cudaGPUDirectRDMAWritesOrderingAllDevices = 200

   ! ---- cudaFlushGPUDirectRDMAWritesScope
   integer(c_int), parameter :: cudaFlushGPUDirectRDMAWritesToOwner = 100
   integer(c_int), parameter :: cudaFlushGPUDirectRDMAWritesToAllDevices = 200

   ! ---- cudaFlushGPUDirectRDMAWritesTarget
   integer(c_int), parameter :: cudaFlushGPUDirectRDMAWritesTargetCurrentDevice = 0

   ! ---- cudaDeviceAttr
   integer(c_int), parameter :: cudaDevAttrMaxThreadsPerBlock = 1
   integer(c_int), parameter :: cudaDevAttrMaxBlockDimX = 2
   integer(c_int), parameter :: cudaDevAttrMaxBlockDimY = 3
   integer(c_int), parameter :: cudaDevAttrMaxBlockDimZ = 4
   integer(c_int), parameter :: cudaDevAttrMaxGridDimX = 5
   integer(c_int), parameter :: cudaDevAttrMaxGridDimY = 6
   integer(c_int), parameter :: cudaDevAttrMaxGridDimZ = 7
   integer(c_int), parameter :: cudaDevAttrMaxSharedMemoryPerBlock = 8
   integer(c_int), parameter :: cudaDevAttrTotalConstantMemory = 9
   integer(c_int), parameter :: cudaDevAttrWarpSize = 10
   integer(c_int), parameter :: cudaDevAttrMaxPitch = 11
   integer(c_int), parameter :: cudaDevAttrMaxRegistersPerBlock = 12
   integer(c_int), parameter :: cudaDevAttrClockRate = 13
   integer(c_int), parameter :: cudaDevAttrTextureAlignment = 14
   integer(c_int), parameter :: cudaDevAttrGpuOverlap = 15
   integer(c_int), parameter :: cudaDevAttrMultiProcessorCount = 16
   integer(c_int), parameter :: cudaDevAttrKernelExecTimeout = 17
   integer(c_int), parameter :: cudaDevAttrIntegrated = 18
   integer(c_int), parameter :: cudaDevAttrCanMapHostMemory = 19
   integer(c_int), parameter :: cudaDevAttrComputeMode = 20
   integer(c_int), parameter :: cudaDevAttrMaxTexture1DWidth = 21
   integer(c_int), parameter :: cudaDevAttrMaxTexture2DWidth = 22
   integer(c_int), parameter :: cudaDevAttrMaxTexture2DHeight = 23
   integer(c_int), parameter :: cudaDevAttrMaxTexture3DWidth = 24
   integer(c_int), parameter :: cudaDevAttrMaxTexture3DHeight = 25
   integer(c_int), parameter :: cudaDevAttrMaxTexture3DDepth = 26
   integer(c_int), parameter :: cudaDevAttrMaxTexture2DLayeredWidth = 27
   integer(c_int), parameter :: cudaDevAttrMaxTexture2DLayeredHeight = 28
   integer(c_int), parameter :: cudaDevAttrMaxTexture2DLayeredLayers = 29
   integer(c_int), parameter :: cudaDevAttrSurfaceAlignment = 30
   integer(c_int), parameter :: cudaDevAttrConcurrentKernels = 31
   integer(c_int), parameter :: cudaDevAttrEccEnabled = 32
   integer(c_int), parameter :: cudaDevAttrPciBusId = 33
   integer(c_int), parameter :: cudaDevAttrPciDeviceId = 34
   integer(c_int), parameter :: cudaDevAttrTccDriver = 35
   integer(c_int), parameter :: cudaDevAttrMemoryClockRate = 36
   integer(c_int), parameter :: cudaDevAttrGlobalMemoryBusWidth = 37
   integer(c_int), parameter :: cudaDevAttrL2CacheSize = 38
   integer(c_int), parameter :: cudaDevAttrMaxThreadsPerMultiProcessor = 39
   integer(c_int), parameter :: cudaDevAttrAsyncEngineCount = 40
   integer(c_int), parameter :: cudaDevAttrUnifiedAddressing = 41
   integer(c_int), parameter :: cudaDevAttrMaxTexture1DLayeredWidth = 42
   integer(c_int), parameter :: cudaDevAttrMaxTexture1DLayeredLayers = 43
   integer(c_int), parameter :: cudaDevAttrMaxTexture2DGatherWidth = 45
   integer(c_int), parameter :: cudaDevAttrMaxTexture2DGatherHeight = 46
   integer(c_int), parameter :: cudaDevAttrMaxTexture3DWidthAlt = 47
   integer(c_int), parameter :: cudaDevAttrMaxTexture3DHeightAlt = 48
   integer(c_int), parameter :: cudaDevAttrMaxTexture3DDepthAlt = 49
   integer(c_int), parameter :: cudaDevAttrPciDomainId = 50
   integer(c_int), parameter :: cudaDevAttrTexturePitchAlignment = 51
   integer(c_int), parameter :: cudaDevAttrMaxTextureCubemapWidth = 52
   integer(c_int), parameter :: cudaDevAttrMaxTextureCubemapLayeredWidth = 53
   integer(c_int), parameter :: cudaDevAttrMaxTextureCubemapLayeredLayers = 54
   integer(c_int), parameter :: cudaDevAttrMaxSurface1DWidth = 55
   integer(c_int), parameter :: cudaDevAttrMaxSurface2DWidth = 56
   integer(c_int), parameter :: cudaDevAttrMaxSurface2DHeight = 57
   integer(c_int), parameter :: cudaDevAttrMaxSurface3DWidth = 58
   integer(c_int), parameter :: cudaDevAttrMaxSurface3DHeight = 59
   integer(c_int), parameter :: cudaDevAttrMaxSurface3DDepth = 60
   integer(c_int), parameter :: cudaDevAttrMaxSurface1DLayeredWidth = 61
   integer(c_int), parameter :: cudaDevAttrMaxSurface1DLayeredLayers = 62
   integer(c_int), parameter :: cudaDevAttrMaxSurface2DLayeredWidth = 63
   integer(c_int), parameter :: cudaDevAttrMaxSurface2DLayeredHeight = 64
   integer(c_int), parameter :: cudaDevAttrMaxSurface2DLayeredLayers = 65
   integer(c_int), parameter :: cudaDevAttrMaxSurfaceCubemapWidth = 66
   integer(c_int), parameter :: cudaDevAttrMaxSurfaceCubemapLayeredWidth = 67
   integer(c_int), parameter :: cudaDevAttrMaxSurfaceCubemapLayeredLayers = 68
   integer(c_int), parameter :: cudaDevAttrMaxTexture1DLinearWidth = 69
   integer(c_int), parameter :: cudaDevAttrMaxTexture2DLinearWidth = 70
   integer(c_int), parameter :: cudaDevAttrMaxTexture2DLinearHeight = 71
   integer(c_int), parameter :: cudaDevAttrMaxTexture2DLinearPitch = 72
   integer(c_int), parameter :: cudaDevAttrMaxTexture2DMipmappedWidth = 73
   integer(c_int), parameter :: cudaDevAttrMaxTexture2DMipmappedHeight = 74
   integer(c_int), parameter :: cudaDevAttrComputeCapabilityMajor = 75
   integer(c_int), parameter :: cudaDevAttrComputeCapabilityMinor = 76
   integer(c_int), parameter :: cudaDevAttrMaxTexture1DMipmappedWidth = 77
   integer(c_int), parameter :: cudaDevAttrStreamPrioritiesSupported = 78
   integer(c_int), parameter :: cudaDevAttrGlobalL1CacheSupported = 79
   integer(c_int), parameter :: cudaDevAttrLocalL1CacheSupported = 80
   integer(c_int), parameter :: cudaDevAttrMaxSharedMemoryPerMultiprocessor = 81
   integer(c_int), parameter :: cudaDevAttrMaxRegistersPerMultiprocessor = 82
   integer(c_int), parameter :: cudaDevAttrManagedMemory = 83
   integer(c_int), parameter :: cudaDevAttrIsMultiGpuBoard = 84
   integer(c_int), parameter :: cudaDevAttrMultiGpuBoardGroupID = 85
   integer(c_int), parameter :: cudaDevAttrHostNativeAtomicSupported = 86
   integer(c_int), parameter :: cudaDevAttrSingleToDoublePrecisionPerfRatio = 87
   integer(c_int), parameter :: cudaDevAttrPageableMemoryAccess = 88
   integer(c_int), parameter :: cudaDevAttrConcurrentManagedAccess = 89
   integer(c_int), parameter :: cudaDevAttrComputePreemptionSupported = 90
   integer(c_int), parameter :: cudaDevAttrCanUseHostPointerForRegisteredMem = 91
   integer(c_int), parameter :: cudaDevAttrReserved92 = 92
   integer(c_int), parameter :: cudaDevAttrReserved93 = 93
   integer(c_int), parameter :: cudaDevAttrReserved94 = 94
   integer(c_int), parameter :: cudaDevAttrCooperativeLaunch = 95
   integer(c_int), parameter :: cudaDevAttrCooperativeMultiDeviceLaunch = 96
   integer(c_int), parameter :: cudaDevAttrMaxSharedMemoryPerBlockOptin = 97
   integer(c_int), parameter :: cudaDevAttrCanFlushRemoteWrites = 98
   integer(c_int), parameter :: cudaDevAttrHostRegisterSupported = 99
   integer(c_int), parameter :: cudaDevAttrPageableMemoryAccessUsesHostPageTables = 100
   integer(c_int), parameter :: cudaDevAttrDirectManagedMemAccessFromHost = 101
   integer(c_int), parameter :: cudaDevAttrMaxBlocksPerMultiprocessor = 106
   integer(c_int), parameter :: cudaDevAttrMaxPersistingL2CacheSize = 108
   integer(c_int), parameter :: cudaDevAttrMaxAccessPolicyWindowSize = 109
   integer(c_int), parameter :: cudaDevAttrReservedSharedMemoryPerBlock = 111
   integer(c_int), parameter :: cudaDevAttrSparseCudaArraySupported = 112
   integer(c_int), parameter :: cudaDevAttrHostRegisterReadOnlySupported = 113
   integer(c_int), parameter :: cudaDevAttrTimelineSemaphoreInteropSupported = 114
   integer(c_int), parameter :: cudaDevAttrMaxTimelineSemaphoreInteropSupported = 114
   integer(c_int), parameter :: cudaDevAttrMemoryPoolsSupported = 115
   integer(c_int), parameter :: cudaDevAttrGPUDirectRDMASupported = 116
   integer(c_int), parameter :: cudaDevAttrGPUDirectRDMAFlushWritesOptions = 117
   integer(c_int), parameter :: cudaDevAttrGPUDirectRDMAWritesOrdering = 118
   integer(c_int), parameter :: cudaDevAttrMemoryPoolSupportedHandleTypes = 119
   integer(c_int), parameter :: cudaDevAttrClusterLaunch = 120
   integer(c_int), parameter :: cudaDevAttrDeferredMappingCudaArraySupported = 121
   integer(c_int), parameter :: cudaDevAttrReserved122 = 122
   integer(c_int), parameter :: cudaDevAttrReserved123 = 123
   integer(c_int), parameter :: cudaDevAttrReserved124 = 124
   integer(c_int), parameter :: cudaDevAttrIpcEventSupport = 125
   integer(c_int), parameter :: cudaDevAttrMemSyncDomainCount = 126
   integer(c_int), parameter :: cudaDevAttrReserved127 = 127
   integer(c_int), parameter :: cudaDevAttrReserved128 = 128
   integer(c_int), parameter :: cudaDevAttrReserved129 = 129
   integer(c_int), parameter :: cudaDevAttrNumaConfig = 130
   integer(c_int), parameter :: cudaDevAttrNumaId = 131
   integer(c_int), parameter :: cudaDevAttrReserved132 = 132
   integer(c_int), parameter :: cudaDevAttrMpsEnabled = 133
   integer(c_int), parameter :: cudaDevAttrHostNumaId = 134
   integer(c_int), parameter :: cudaDevAttrD3D12CigSupported = 135
   integer(c_int), parameter :: cudaDevAttrVulkanCigSupported = 138
   integer(c_int), parameter :: cudaDevAttrGpuPciDeviceId = 139
   integer(c_int), parameter :: cudaDevAttrGpuPciSubsystemId = 140
   integer(c_int), parameter :: cudaDevAttrReserved141 = 141
   integer(c_int), parameter :: cudaDevAttrHostNumaMemoryPoolsSupported = 142
   integer(c_int), parameter :: cudaDevAttrHostNumaMultinodeIpcSupported = 143
   integer(c_int), parameter :: cudaDevAttrMax = 144

   ! ---- cudaMemPoolAttr
   integer(c_int), parameter :: cudaMemPoolReuseFollowEventDependencies = 1
   integer(c_int), parameter :: cudaMemPoolReuseAllowOpportunistic = 2
   integer(c_int), parameter :: cudaMemPoolReuseAllowInternalDependencies = 3
   integer(c_int), parameter :: cudaMemPoolAttrReleaseThreshold = 4
   integer(c_int), parameter :: cudaMemPoolAttrReservedMemCurrent = 5
   integer(c_int), parameter :: cudaMemPoolAttrReservedMemHigh = 6
   integer(c_int), parameter :: cudaMemPoolAttrUsedMemCurrent = 7
   integer(c_int), parameter :: cudaMemPoolAttrUsedMemHigh = 8

   ! ---- cudaMemLocationType
   integer(c_int), parameter :: cudaMemLocationTypeInvalid = 0
   integer(c_int), parameter :: cudaMemLocationTypeDevice = 1
   integer(c_int), parameter :: cudaMemLocationTypeHost = 2
   integer(c_int), parameter :: cudaMemLocationTypeHostNuma = 3
   integer(c_int), parameter :: cudaMemLocationTypeHostNumaCurrent = 4

   ! ---- cudaMemAccessFlags
   integer(c_int), parameter :: cudaMemAccessFlagsProtNone = 0
   integer(c_int), parameter :: cudaMemAccessFlagsProtRead = 1
   integer(c_int), parameter :: cudaMemAccessFlagsProtReadWrite = 3

   ! ---- cudaMemAllocationType
   integer(c_int), parameter :: cudaMemAllocationTypeInvalid = 0
   integer(c_int), parameter :: cudaMemAllocationTypePinned = 1
   integer(c_int), parameter :: cudaMemAllocationTypeMax = 2147483647

   ! ---- cudaMemAllocationHandleType
   integer(c_int), parameter :: cudaMemHandleTypeNone = 0
   integer(c_int), parameter :: cudaMemHandleTypePosixFileDescriptor = 1
   integer(c_int), parameter :: cudaMemHandleTypeWin32 = 2
   integer(c_int), parameter :: cudaMemHandleTypeWin32Kmt = 4
   integer(c_int), parameter :: cudaMemHandleTypeFabric = 8

   ! ---- cudaGraphMemAttributeType
   integer(c_int), parameter :: cudaGraphMemAttrUsedMemCurrent = 0
   integer(c_int), parameter :: cudaGraphMemAttrUsedMemHigh = 1
   integer(c_int), parameter :: cudaGraphMemAttrReservedMemCurrent = 2
   integer(c_int), parameter :: cudaGraphMemAttrReservedMemHigh = 3

   ! ---- cudaMemcpyFlags
   integer(c_int), parameter :: cudaMemcpyFlagDefault = 0
   integer(c_int), parameter :: cudaMemcpyFlagPreferOverlapWithCompute = 1

   ! ---- cudaMemcpySrcAccessOrder
   integer(c_int), parameter :: cudaMemcpySrcAccessOrderInvalid = 0
   integer(c_int), parameter :: cudaMemcpySrcAccessOrderStream = 1
   integer(c_int), parameter :: cudaMemcpySrcAccessOrderDuringApiCall = 2
   integer(c_int), parameter :: cudaMemcpySrcAccessOrderAny = 3
   integer(c_int), parameter :: cudaMemcpySrcAccessOrderMax = 2147483647

   ! ---- cudaMemcpy3DOperandType
   integer(c_int), parameter :: cudaMemcpyOperandTypePointer = 1
   integer(c_int), parameter :: cudaMemcpyOperandTypeArray = 2
   integer(c_int), parameter :: cudaMemcpyOperandTypeMax = 2147483647

   ! ---- cudaDeviceP2PAttr
   integer(c_int), parameter :: cudaDevP2PAttrPerformanceRank = 1
   integer(c_int), parameter :: cudaDevP2PAttrAccessSupported = 2
   integer(c_int), parameter :: cudaDevP2PAttrNativeAtomicSupported = 3
   integer(c_int), parameter :: cudaDevP2PAttrCudaArrayAccessSupported = 4

   ! ---- cudaExternalMemoryHandleType
   integer(c_int), parameter :: cudaExternalMemoryHandleTypeOpaqueFd = 1
   integer(c_int), parameter :: cudaExternalMemoryHandleTypeOpaqueWin32 = 2
   integer(c_int), parameter :: cudaExternalMemoryHandleTypeOpaqueWin32Kmt = 3
   integer(c_int), parameter :: cudaExternalMemoryHandleTypeD3D12Heap = 4
   integer(c_int), parameter :: cudaExternalMemoryHandleTypeD3D12Resource = 5
   integer(c_int), parameter :: cudaExternalMemoryHandleTypeD3D11Resource = 6
   integer(c_int), parameter :: cudaExternalMemoryHandleTypeD3D11ResourceKmt = 7
   integer(c_int), parameter :: cudaExternalMemoryHandleTypeNvSciBuf = 8

   ! ---- cudaExternalSemaphoreHandleType
   integer(c_int), parameter :: cudaExternalSemaphoreHandleTypeOpaqueFd = 1
   integer(c_int), parameter :: cudaExternalSemaphoreHandleTypeOpaqueWin32 = 2
   integer(c_int), parameter :: cudaExternalSemaphoreHandleTypeOpaqueWin32Kmt = 3
   integer(c_int), parameter :: cudaExternalSemaphoreHandleTypeD3D12Fence = 4
   integer(c_int), parameter :: cudaExternalSemaphoreHandleTypeD3D11Fence = 5
   integer(c_int), parameter :: cudaExternalSemaphoreHandleTypeNvSciSync = 6
   integer(c_int), parameter :: cudaExternalSemaphoreHandleTypeKeyedMutex = 7
   integer(c_int), parameter :: cudaExternalSemaphoreHandleTypeKeyedMutexKmt = 8
   integer(c_int), parameter :: cudaExternalSemaphoreHandleTypeTimelineSemaphoreFd = 9
   integer(c_int), parameter :: cudaExternalSemaphoreHandleTypeTimelineSemaphoreWin32 = 10

   ! ---- cudaJitOption
   integer(c_int), parameter :: cudaJitMaxRegisters = 0
   integer(c_int), parameter :: cudaJitThreadsPerBlock = 1
   integer(c_int), parameter :: cudaJitWallTime = 2
   integer(c_int), parameter :: cudaJitInfoLogBuffer = 3
   integer(c_int), parameter :: cudaJitInfoLogBufferSizeBytes = 4
   integer(c_int), parameter :: cudaJitErrorLogBuffer = 5
   integer(c_int), parameter :: cudaJitErrorLogBufferSizeBytes = 6
   integer(c_int), parameter :: cudaJitOptimizationLevel = 7
   integer(c_int), parameter :: cudaJitFallbackStrategy = 10
   integer(c_int), parameter :: cudaJitGenerateDebugInfo = 11
   integer(c_int), parameter :: cudaJitLogVerbose = 12
   integer(c_int), parameter :: cudaJitGenerateLineInfo = 13
   integer(c_int), parameter :: cudaJitCacheMode = 14
   integer(c_int), parameter :: cudaJitPositionIndependentCode = 30
   integer(c_int), parameter :: cudaJitMinCtaPerSm = 31
   integer(c_int), parameter :: cudaJitMaxThreadsPerBlock = 32
   integer(c_int), parameter :: cudaJitOverrideDirectiveValues = 33

   ! ---- cudaLibraryOption
   integer(c_int), parameter :: cudaLibraryHostUniversalFunctionAndDataTable = 0
   integer(c_int), parameter :: cudaLibraryBinaryIsPreserved = 1

   ! ---- cudaJit_CacheMode
   integer(c_int), parameter :: cudaJitCacheOptionNone = 0
   integer(c_int), parameter :: cudaJitCacheOptionCG = 1
   integer(c_int), parameter :: cudaJitCacheOptionCA = 2

   ! ---- cudaJit_Fallback
   integer(c_int), parameter :: cudaPreferPtx = 0
   integer(c_int), parameter :: cudaPreferBinary = 1

   ! ---- cudaCGScope
   integer(c_int), parameter :: cudaCGScopeInvalid = 0
   integer(c_int), parameter :: cudaCGScopeGrid = 1
   integer(c_int), parameter :: cudaCGScopeMultiGrid = 2

   ! ---- cudaGraphConditionalHandleFlags
   integer(c_int), parameter :: cudaGraphCondAssignDefault = 1

   ! ---- cudaGraphConditionalNodeType
   integer(c_int), parameter :: cudaGraphCondTypeIf = 0
   integer(c_int), parameter :: cudaGraphCondTypeWhile = 1
   integer(c_int), parameter :: cudaGraphCondTypeSwitch = 2

   ! ---- cudaGraphNodeType
   integer(c_int), parameter :: cudaGraphNodeTypeKernel = 0
   integer(c_int), parameter :: cudaGraphNodeTypeMemcpy = 1
   integer(c_int), parameter :: cudaGraphNodeTypeMemset = 2
   integer(c_int), parameter :: cudaGraphNodeTypeHost = 3
   integer(c_int), parameter :: cudaGraphNodeTypeGraph = 4
   integer(c_int), parameter :: cudaGraphNodeTypeEmpty = 5
   integer(c_int), parameter :: cudaGraphNodeTypeWaitEvent = 6
   integer(c_int), parameter :: cudaGraphNodeTypeEventRecord = 7
   integer(c_int), parameter :: cudaGraphNodeTypeExtSemaphoreSignal = 8
   integer(c_int), parameter :: cudaGraphNodeTypeExtSemaphoreWait = 9
   integer(c_int), parameter :: cudaGraphNodeTypeMemAlloc = 10
   integer(c_int), parameter :: cudaGraphNodeTypeMemFree = 11
   integer(c_int), parameter :: cudaGraphNodeTypeConditional = 13
   integer(c_int), parameter :: cudaGraphNodeTypeCount = 14

   ! ---- cudaGraphChildGraphNodeOwnership
   integer(c_int), parameter :: cudaGraphChildGraphOwnershipClone = 0
   integer(c_int), parameter :: cudaGraphChildGraphOwnershipMove = 1

   ! ---- cudaGraphDependencyType
   integer(c_int), parameter :: cudaGraphDependencyTypeDefault = 0
   integer(c_int), parameter :: cudaGraphDependencyTypeProgrammatic = 1

   ! ---- cudaGraphExecUpdateResult
   integer(c_int), parameter :: cudaGraphExecUpdateSuccess = 0
   integer(c_int), parameter :: cudaGraphExecUpdateError = 1
   integer(c_int), parameter :: cudaGraphExecUpdateErrorTopologyChanged = 2
   integer(c_int), parameter :: cudaGraphExecUpdateErrorNodeTypeChanged = 3
   integer(c_int), parameter :: cudaGraphExecUpdateErrorFunctionChanged = 4
   integer(c_int), parameter :: cudaGraphExecUpdateErrorParametersChanged = 5
   integer(c_int), parameter :: cudaGraphExecUpdateErrorNotSupported = 6
   integer(c_int), parameter :: cudaGraphExecUpdateErrorUnsupportedFunctionChange = 7
   integer(c_int), parameter :: cudaGraphExecUpdateErrorAttributesChanged = 8

   ! ---- cudaGraphInstantiateResult
   integer(c_int), parameter :: cudaGraphInstantiateSuccess = 0
   integer(c_int), parameter :: cudaGraphInstantiateError = 1
   integer(c_int), parameter :: cudaGraphInstantiateInvalidStructure = 2
   integer(c_int), parameter :: cudaGraphInstantiateNodeOperationNotSupported = 3
   integer(c_int), parameter :: cudaGraphInstantiateMultipleDevicesNotSupported = 4
   integer(c_int), parameter :: cudaGraphInstantiateConditionalHandleUnused = 5

   ! ---- cudaGraphKernelNodeField
   integer(c_int), parameter :: cudaGraphKernelNodeFieldInvalid = 0
   integer(c_int), parameter :: cudaGraphKernelNodeFieldGridDim = 1
   integer(c_int), parameter :: cudaGraphKernelNodeFieldParam = 2
   integer(c_int), parameter :: cudaGraphKernelNodeFieldEnabled = 3

   ! ---- cudaGetDriverEntryPointFlags
   integer(c_int), parameter :: cudaEnableDefault = 0
   integer(c_int), parameter :: cudaEnableLegacyStream = 1
   integer(c_int), parameter :: cudaEnablePerThreadDefaultStream = 2

   ! ---- cudaDriverEntryPointQueryResult
   integer(c_int), parameter :: cudaDriverEntryPointSuccess = 0
   integer(c_int), parameter :: cudaDriverEntryPointSymbolNotFound = 1
   integer(c_int), parameter :: cudaDriverEntryPointVersionNotSufficent = 2

   ! ---- cudaGraphDebugDotFlags
   integer(c_int), parameter :: cudaGraphDebugDotFlagsVerbose = 1
   integer(c_int), parameter :: cudaGraphDebugDotFlagsKernelNodeParams = 4
   integer(c_int), parameter :: cudaGraphDebugDotFlagsMemcpyNodeParams = 8
   integer(c_int), parameter :: cudaGraphDebugDotFlagsMemsetNodeParams = 16
   integer(c_int), parameter :: cudaGraphDebugDotFlagsHostNodeParams = 32
   integer(c_int), parameter :: cudaGraphDebugDotFlagsEventNodeParams = 64
   integer(c_int), parameter :: cudaGraphDebugDotFlagsExtSemasSignalNodeParams = 128
   integer(c_int), parameter :: cudaGraphDebugDotFlagsExtSemasWaitNodeParams = 256
   integer(c_int), parameter :: cudaGraphDebugDotFlagsKernelNodeAttributes = 512
   integer(c_int), parameter :: cudaGraphDebugDotFlagsHandles = 1024
   integer(c_int), parameter :: cudaGraphDebugDotFlagsConditionalNodeParams = 32768

   ! ---- cudaGraphInstantiateFlags
   integer(c_int), parameter :: cudaGraphInstantiateFlagAutoFreeOnLaunch = 1
   integer(c_int), parameter :: cudaGraphInstantiateFlagUpload = 2
   integer(c_int), parameter :: cudaGraphInstantiateFlagDeviceLaunch = 4
   integer(c_int), parameter :: cudaGraphInstantiateFlagUseNodePriority = 8

   ! ---- cudaLaunchMemSyncDomain
   integer(c_int), parameter :: cudaLaunchMemSyncDomainDefault = 0
   integer(c_int), parameter :: cudaLaunchMemSyncDomainRemote = 1

   ! ---- cudaLaunchAttributeID
   integer(c_int), parameter :: cudaLaunchAttributeIgnore = 0
   integer(c_int), parameter :: cudaLaunchAttributeAccessPolicyWindow = 1
   integer(c_int), parameter :: cudaLaunchAttributeCooperative = 2
   integer(c_int), parameter :: cudaLaunchAttributeSynchronizationPolicy = 3
   integer(c_int), parameter :: cudaLaunchAttributeClusterDimension = 4
   integer(c_int), parameter :: cudaLaunchAttributeClusterSchedulingPolicyPreference = 5
   integer(c_int), parameter :: cudaLaunchAttributeProgrammaticStreamSerialization = 6
   integer(c_int), parameter :: cudaLaunchAttributeProgrammaticEvent = 7
   integer(c_int), parameter :: cudaLaunchAttributePriority = 8
   integer(c_int), parameter :: cudaLaunchAttributeMemSyncDomainMap = 9
   integer(c_int), parameter :: cudaLaunchAttributeMemSyncDomain = 10
   integer(c_int), parameter :: cudaLaunchAttributePreferredClusterDimension = 11
   integer(c_int), parameter :: cudaLaunchAttributeLaunchCompletionEvent = 12
   integer(c_int), parameter :: cudaLaunchAttributeDeviceUpdatableKernelNode = 13
   integer(c_int), parameter :: cudaLaunchAttributePreferredSharedMemoryCarveout = 14

   ! ---- cudaDeviceNumaConfig
   integer(c_int), parameter :: cudaDeviceNumaConfigNone = 0
   integer(c_int), parameter :: cudaDeviceNumaConfigNumaNode = 1

   ! ---- cudaAsyncNotificationType
   integer(c_int), parameter :: cudaAsyncNotificationTypeOverBudget = 1

   ! ---- cudaSurfaceBoundaryMode
   integer(c_int), parameter :: cudaBoundaryModeZero = 0
   integer(c_int), parameter :: cudaBoundaryModeClamp = 1
   integer(c_int), parameter :: cudaBoundaryModeTrap = 2

   ! ---- cudaSurfaceFormatMode
   integer(c_int), parameter :: cudaFormatModeForced = 0
   integer(c_int), parameter :: cudaFormatModeAuto = 1

   ! ---- cudaTextureAddressMode
   integer(c_int), parameter :: cudaAddressModeWrap = 0
   integer(c_int), parameter :: cudaAddressModeClamp = 1
   integer(c_int), parameter :: cudaAddressModeMirror = 2
   integer(c_int), parameter :: cudaAddressModeBorder = 3

   ! ---- cudaTextureFilterMode
   integer(c_int), parameter :: cudaFilterModePoint = 0
   integer(c_int), parameter :: cudaFilterModeLinear = 1

   ! ---- cudaTextureReadMode
   integer(c_int), parameter :: cudaReadModeElementType = 0
   integer(c_int), parameter :: cudaReadModeNormalizedFloat = 1

   ! ---- cudaError

   ! ---- cudaChannelFormatKind

   ! ---- cudaMemoryType

   ! ---- cudaMemcpyKind

   ! ---- cudaAccessProperty

   ! ---- cudaStreamCaptureStatus

   ! ---- cudaStreamCaptureMode

   ! ---- cudaSynchronizationPolicy

   ! ---- cudaClusterSchedulingPolicy

   ! ---- cudaStreamUpdateCaptureDependenciesFlags

   ! ---- cudaUserObjectFlags

   ! ---- cudaUserObjectRetainFlags

   ! ---- cudaGraphicsRegisterFlags

   ! ---- cudaGraphicsMapFlags

   ! ---- cudaGraphicsCubeFace

   ! ---- cudaResourceType

   ! ---- cudaResourceViewFormat

   ! ---- cudaFuncAttribute

   ! ---- cudaFuncCache

   ! ---- cudaSharedMemConfig

   ! ---- cudaSharedCarveout

   ! ---- cudaComputeMode

   ! ---- cudaLimit

   ! ---- cudaMemoryAdvise

   ! ---- cudaMemRangeAttribute

   ! ---- cudaFlushGPUDirectRDMAWritesOptions

   ! ---- cudaGPUDirectRDMAWritesOrdering

   ! ---- cudaFlushGPUDirectRDMAWritesScope

   ! ---- cudaFlushGPUDirectRDMAWritesTarget

   ! ---- cudaDeviceAttr

   ! ---- cudaMemPoolAttr

   ! ---- cudaMemLocationType

   ! ---- cudaMemAccessFlags

   ! ---- cudaMemAllocationType

   ! ---- cudaMemAllocationHandleType

   ! ---- cudaGraphMemAttributeType

   ! ---- cudaMemcpyFlags

   ! ---- cudaMemcpySrcAccessOrder

   ! ---- cudaMemcpy3DOperandType

   ! ---- cudaDeviceP2PAttr

   ! ---- cudaExternalMemoryHandleType

   ! ---- cudaExternalSemaphoreHandleType

   ! ---- cudaJitOption

   ! ---- cudaLibraryOption

   ! ---- cudaJit_CacheMode

   ! ---- cudaJit_Fallback

   ! ---- cudaCGScope

   ! ---- cudaGraphConditionalHandleFlags

   ! ---- cudaGraphConditionalNodeType

   ! ---- cudaGraphNodeType

   ! ---- cudaGraphChildGraphNodeOwnership

   ! ---- cudaGraphDependencyType

   ! ---- cudaGraphExecUpdateResult

   ! ---- cudaGraphInstantiateResult

   ! ---- cudaGraphKernelNodeField

   ! ---- cudaGetDriverEntryPointFlags

   ! ---- cudaDriverEntryPointQueryResult

   ! ---- cudaGraphDebugDotFlags

   ! ---- cudaGraphInstantiateFlags

   ! ---- cudaLaunchMemSyncDomain

   ! ---- cudaLaunchAttributeID

   ! ---- cudaDeviceNumaConfig

   ! ---- cudaAsyncNotificationType

   ! ---- cudaDataType
   integer(c_int), parameter :: CUDA_R_16F = 2
   integer(c_int), parameter :: CUDA_C_16F = 6
   integer(c_int), parameter :: CUDA_R_16BF = 14
   integer(c_int), parameter :: CUDA_C_16BF = 15
   integer(c_int), parameter :: CUDA_R_32F = 0
   integer(c_int), parameter :: CUDA_C_32F = 4
   integer(c_int), parameter :: CUDA_R_64F = 1
   integer(c_int), parameter :: CUDA_C_64F = 5
   integer(c_int), parameter :: CUDA_R_4I = 16
   integer(c_int), parameter :: CUDA_C_4I = 17
   integer(c_int), parameter :: CUDA_R_4U = 18
   integer(c_int), parameter :: CUDA_C_4U = 19
   integer(c_int), parameter :: CUDA_R_8I = 3
   integer(c_int), parameter :: CUDA_C_8I = 7
   integer(c_int), parameter :: CUDA_R_8U = 8
   integer(c_int), parameter :: CUDA_C_8U = 9
   integer(c_int), parameter :: CUDA_R_16I = 20
   integer(c_int), parameter :: CUDA_C_16I = 21
   integer(c_int), parameter :: CUDA_R_16U = 22
   integer(c_int), parameter :: CUDA_C_16U = 23
   integer(c_int), parameter :: CUDA_R_32I = 10
   integer(c_int), parameter :: CUDA_C_32I = 11
   integer(c_int), parameter :: CUDA_R_32U = 12
   integer(c_int), parameter :: CUDA_C_32U = 13
   integer(c_int), parameter :: CUDA_R_64I = 24
   integer(c_int), parameter :: CUDA_C_64I = 25
   integer(c_int), parameter :: CUDA_R_64U = 26
   integer(c_int), parameter :: CUDA_C_64U = 27
   integer(c_int), parameter :: CUDA_R_8F_E4M3 = 28
   integer(c_int), parameter :: CUDA_R_8F_E5M2 = 29
   integer(c_int), parameter :: CUDA_R_8F_UE8M0 = 30
   integer(c_int), parameter :: CUDA_R_6F_E2M3 = 31
   integer(c_int), parameter :: CUDA_R_6F_E3M2 = 32
   integer(c_int), parameter :: CUDA_R_4F_E2M1 = 33

   ! ---- libraryPropertyType
   integer(c_int), parameter :: MAJOR_VERSION = 0
   integer(c_int), parameter :: MINOR_VERSION = 1
   integer(c_int), parameter :: PATCH_LEVEL = 2

   ! ---- cudaError

   ! ---- cudaChannelFormatKind

   ! ---- cudaMemoryType

   ! ---- cudaMemcpyKind

   ! ---- cudaAccessProperty

   ! ---- cudaStreamCaptureStatus

   ! ---- cudaStreamCaptureMode

   ! ---- cudaSynchronizationPolicy

   ! ---- cudaClusterSchedulingPolicy

   ! ---- cudaStreamUpdateCaptureDependenciesFlags

   ! ---- cudaUserObjectFlags

   ! ---- cudaUserObjectRetainFlags

   ! ---- cudaGraphicsRegisterFlags

   ! ---- cudaGraphicsMapFlags

   ! ---- cudaGraphicsCubeFace

   ! ---- cudaResourceType

   ! ---- cudaResourceViewFormat

   ! ---- cudaFuncAttribute

   ! ---- cudaFuncCache

   ! ---- cudaSharedMemConfig

   ! ---- cudaSharedCarveout

   ! ---- cudaComputeMode

   ! ---- cudaLimit

   ! ---- cudaMemoryAdvise

   ! ---- cudaMemRangeAttribute

   ! ---- cudaFlushGPUDirectRDMAWritesOptions

   ! ---- cudaGPUDirectRDMAWritesOrdering

   ! ---- cudaFlushGPUDirectRDMAWritesScope

   ! ---- cudaFlushGPUDirectRDMAWritesTarget

   ! ---- cudaDeviceAttr

   ! ---- cudaMemPoolAttr

   ! ---- cudaMemLocationType

   ! ---- cudaMemAccessFlags

   ! ---- cudaMemAllocationType

   ! ---- cudaMemAllocationHandleType

   ! ---- cudaGraphMemAttributeType

   ! ---- cudaMemcpyFlags

   ! ---- cudaMemcpySrcAccessOrder

   ! ---- cudaMemcpy3DOperandType

   ! ---- cudaDeviceP2PAttr

   ! ---- cudaExternalMemoryHandleType

   ! ---- cudaExternalSemaphoreHandleType

   ! ---- cudaJitOption

   ! ---- cudaLibraryOption

   ! ---- cudaJit_CacheMode

   ! ---- cudaJit_Fallback

   ! ---- cudaCGScope

   ! ---- cudaGraphConditionalHandleFlags

   ! ---- cudaGraphConditionalNodeType

   ! ---- cudaGraphNodeType

   ! ---- cudaGraphChildGraphNodeOwnership

   ! ---- cudaGraphDependencyType

   ! ---- cudaGraphExecUpdateResult

   ! ---- cudaGraphInstantiateResult

   ! ---- cudaGraphKernelNodeField

   ! ---- cudaGetDriverEntryPointFlags

   ! ---- cudaDriverEntryPointQueryResult

   ! ---- cudaGraphDebugDotFlags

   ! ---- cudaGraphInstantiateFlags

   ! ---- cudaLaunchMemSyncDomain

   ! ---- cudaLaunchAttributeID

   ! ---- cudaDeviceNumaConfig

   ! ---- cudaAsyncNotificationType

   ! ---- cudaTextureAddressMode

   ! ---- cudaTextureFilterMode

   ! ---- cudaTextureReadMode

   ! ---- cudaError

   ! ---- cudaChannelFormatKind

   ! ---- cudaMemoryType

   ! ---- cudaMemcpyKind

   ! ---- cudaAccessProperty

   ! ---- cudaStreamCaptureStatus

   ! ---- cudaStreamCaptureMode

   ! ---- cudaSynchronizationPolicy

   ! ---- cudaClusterSchedulingPolicy

   ! ---- cudaStreamUpdateCaptureDependenciesFlags

   ! ---- cudaUserObjectFlags

   ! ---- cudaUserObjectRetainFlags

   ! ---- cudaGraphicsRegisterFlags

   ! ---- cudaGraphicsMapFlags

   ! ---- cudaGraphicsCubeFace

   ! ---- cudaResourceType

   ! ---- cudaResourceViewFormat

   ! ---- cudaFuncAttribute

   ! ---- cudaFuncCache

   ! ---- cudaSharedMemConfig

   ! ---- cudaSharedCarveout

   ! ---- cudaComputeMode

   ! ---- cudaLimit

   ! ---- cudaMemoryAdvise

   ! ---- cudaMemRangeAttribute

   ! ---- cudaFlushGPUDirectRDMAWritesOptions

   ! ---- cudaGPUDirectRDMAWritesOrdering

   ! ---- cudaFlushGPUDirectRDMAWritesScope

   ! ---- cudaFlushGPUDirectRDMAWritesTarget

   ! ---- cudaDeviceAttr

   ! ---- cudaMemPoolAttr

   ! ---- cudaMemLocationType

   ! ---- cudaMemAccessFlags

   ! ---- cudaMemAllocationType

   ! ---- cudaMemAllocationHandleType

   ! ---- cudaGraphMemAttributeType

   ! ---- cudaMemcpyFlags

   ! ---- cudaMemcpySrcAccessOrder

   ! ---- cudaMemcpy3DOperandType

   ! ---- cudaDeviceP2PAttr

   ! ---- cudaExternalMemoryHandleType

   ! ---- cudaExternalSemaphoreHandleType

   ! ---- cudaJitOption

   ! ---- cudaLibraryOption

   ! ---- cudaJit_CacheMode

   ! ---- cudaJit_Fallback

   ! ---- cudaCGScope

   ! ---- cudaGraphConditionalHandleFlags

   ! ---- cudaGraphConditionalNodeType

   ! ---- cudaGraphNodeType

   ! ---- cudaGraphChildGraphNodeOwnership

   ! ---- cudaGraphDependencyType

   ! ---- cudaGraphExecUpdateResult

   ! ---- cudaGraphInstantiateResult

   ! ---- cudaGraphKernelNodeField

   ! ---- cudaGetDriverEntryPointFlags

   ! ---- cudaDriverEntryPointQueryResult

   ! ---- cudaGraphDebugDotFlags

   ! ---- cudaGraphInstantiateFlags

   ! ---- cudaLaunchMemSyncDomain

   ! ---- cudaLaunchAttributeID

   ! ---- cudaDeviceNumaConfig

   ! ---- cudaAsyncNotificationType

   ! ---- cudaSurfaceBoundaryMode

   ! ---- cudaSurfaceFormatMode

   ! ---- cudaError

   ! ---- cudaChannelFormatKind

   ! ---- cudaMemoryType

   ! ---- cudaMemcpyKind

   ! ---- cudaAccessProperty

   ! ---- cudaStreamCaptureStatus

   ! ---- cudaStreamCaptureMode

   ! ---- cudaSynchronizationPolicy

   ! ---- cudaClusterSchedulingPolicy

   ! ---- cudaStreamUpdateCaptureDependenciesFlags

   ! ---- cudaUserObjectFlags

   ! ---- cudaUserObjectRetainFlags

   ! ---- cudaGraphicsRegisterFlags

   ! ---- cudaGraphicsMapFlags

   ! ---- cudaGraphicsCubeFace

   ! ---- cudaResourceType

   ! ---- cudaResourceViewFormat

   ! ---- cudaFuncAttribute

   ! ---- cudaFuncCache

   ! ---- cudaSharedMemConfig

   ! ---- cudaSharedCarveout

   ! ---- cudaComputeMode

   ! ---- cudaLimit

   ! ---- cudaMemoryAdvise

   ! ---- cudaMemRangeAttribute

   ! ---- cudaFlushGPUDirectRDMAWritesOptions

   ! ---- cudaGPUDirectRDMAWritesOrdering

   ! ---- cudaFlushGPUDirectRDMAWritesScope

   ! ---- cudaFlushGPUDirectRDMAWritesTarget

   ! ---- cudaDeviceAttr

   ! ---- cudaMemPoolAttr

   ! ---- cudaMemLocationType

   ! ---- cudaMemAccessFlags

   ! ---- cudaMemAllocationType

   ! ---- cudaMemAllocationHandleType

   ! ---- cudaGraphMemAttributeType

   ! ---- cudaMemcpyFlags

   ! ---- cudaMemcpySrcAccessOrder

   ! ---- cudaMemcpy3DOperandType

   ! ---- cudaDeviceP2PAttr

   ! ---- cudaExternalMemoryHandleType

   ! ---- cudaExternalSemaphoreHandleType

   ! ---- cudaJitOption

   ! ---- cudaLibraryOption

   ! ---- cudaJit_CacheMode

   ! ---- cudaJit_Fallback

   ! ---- cudaCGScope

   ! ---- cudaGraphConditionalHandleFlags

   ! ---- cudaGraphConditionalNodeType

   ! ---- cudaGraphNodeType

   ! ---- cudaGraphChildGraphNodeOwnership

   ! ---- cudaGraphDependencyType

   ! ---- cudaGraphExecUpdateResult

   ! ---- cudaGraphInstantiateResult

   ! ---- cudaGraphKernelNodeField

   ! ---- cudaGetDriverEntryPointFlags

   ! ---- cudaDriverEntryPointQueryResult

   ! ---- cudaGraphDebugDotFlags

   ! ---- cudaGraphInstantiateFlags

   ! ---- cudaLaunchMemSyncDomain

   ! ---- cudaLaunchAttributeID

   ! ---- cudaDeviceNumaConfig

   ! ---- cudaAsyncNotificationType

   ! ======================================================================
   !  Flag constants defined as C macros rather than enumerators
   !  (stream/event/host-alloc/device-schedule flags, ...)
   ! ======================================================================
   integer(c_int), parameter :: cudaArrayColorAttachment = 32
   integer(c_int), parameter :: cudaArrayCubemap = 4
   integer(c_int), parameter :: cudaArrayDefault = 0
   integer(c_int), parameter :: cudaArrayDeferredMapping = 128
   integer(c_int), parameter :: cudaArrayLayered = 1
   integer(c_int), parameter :: cudaArraySparse = 64
   integer(c_int), parameter :: cudaArraySparsePropertiesSingleMipTail = 1
   integer(c_int), parameter :: cudaArraySurfaceLoadStore = 2
   integer(c_int), parameter :: cudaArrayTextureGather = 8
   integer(c_int), parameter :: cudaCooperativeLaunchMultiDeviceNoPostSync = 2
   integer(c_int), parameter :: cudaCooperativeLaunchMultiDeviceNoPreSync = 1
   integer(c_int), parameter :: cudaDeviceBlockingSync = 4
   integer(c_int), parameter :: cudaDeviceLmemResizeToMax = 16
   integer(c_int), parameter :: cudaDeviceMapHost = 8
   integer(c_int), parameter :: cudaDeviceMask = 255
   integer(c_int), parameter :: cudaDeviceScheduleAuto = 0
   integer(c_int), parameter :: cudaDeviceScheduleBlockingSync = 4
   integer(c_int), parameter :: cudaDeviceScheduleMask = 7
   integer(c_int), parameter :: cudaDeviceScheduleSpin = 1
   integer(c_int), parameter :: cudaDeviceScheduleYield = 2
   integer(c_int), parameter :: cudaDeviceSyncMemops = 128
   integer(c_int), parameter :: cudaEventBlockingSync = 1
   integer(c_int), parameter :: cudaEventDefault = 0
   integer(c_int), parameter :: cudaEventDisableTiming = 2
   integer(c_int), parameter :: cudaEventInterprocess = 4
   integer(c_int), parameter :: cudaEventRecordDefault = 0
   integer(c_int), parameter :: cudaEventRecordExternal = 1
   integer(c_int), parameter :: cudaEventWaitDefault = 0
   integer(c_int), parameter :: cudaEventWaitExternal = 1
   integer(c_int), parameter :: cudaExternalMemoryDedicated = 1
   integer(c_int), parameter :: cudaExternalSemaphoreSignalSkipNvSciBufMemSync = 1
   integer(c_int), parameter :: cudaExternalSemaphoreWaitSkipNvSciBufMemSync = 2
   integer(c_int), parameter :: cudaGraphKernelNodePortDefault = 0
   integer(c_int), parameter :: cudaGraphKernelNodePortLaunchCompletion = 2
   integer(c_int), parameter :: cudaGraphKernelNodePortProgrammatic = 1
   integer(c_int), parameter :: cudaHostAllocDefault = 0
   integer(c_int), parameter :: cudaHostAllocMapped = 2
   integer(c_int), parameter :: cudaHostAllocPortable = 1
   integer(c_int), parameter :: cudaHostAllocWriteCombined = 4
   integer(c_int), parameter :: cudaHostRegisterDefault = 0
   integer(c_int), parameter :: cudaHostRegisterIoMemory = 4
   integer(c_int), parameter :: cudaHostRegisterMapped = 2
   integer(c_int), parameter :: cudaHostRegisterPortable = 1
   integer(c_int), parameter :: cudaHostRegisterReadOnly = 8
   integer(c_int), parameter :: cudaInitDeviceFlagsAreValid = 1
   integer(c_int), parameter :: cudaIpcMemLazyEnablePeerAccess = 1
   integer(c_int), parameter :: cudaMemAttachGlobal = 1
   integer(c_int), parameter :: cudaMemAttachHost = 2
   integer(c_int), parameter :: cudaMemAttachSingle = 4
   integer(c_int), parameter :: cudaMemPoolCreateUsageHwDecompress = 2
   integer(c_int), parameter :: cudaNvSciSyncAttrSignal = 1
   integer(c_int), parameter :: cudaNvSciSyncAttrWait = 2
   integer(c_int), parameter :: cudaOccupancyDefault = 0
   integer(c_int), parameter :: cudaOccupancyDisableCachingOverride = 1
   integer(c_int), parameter :: cudaPeerAccessDefault = 0
   integer(c_int), parameter :: cudaStreamDefault = 0
   integer(c_int), parameter :: cudaStreamNonBlocking = 1
   integer(c_int), parameter :: cudaSurfaceType1D = 1
   integer(c_int), parameter :: cudaSurfaceType1DLayered = 241
   integer(c_int), parameter :: cudaSurfaceType2D = 2
   integer(c_int), parameter :: cudaSurfaceType2DLayered = 242
   integer(c_int), parameter :: cudaSurfaceType3D = 3
   integer(c_int), parameter :: cudaSurfaceTypeCubemap = 12
   integer(c_int), parameter :: cudaSurfaceTypeCubemapLayered = 252
   integer(c_int), parameter :: cudaTextureType1D = 1
   integer(c_int), parameter :: cudaTextureType1DLayered = 241
   integer(c_int), parameter :: cudaTextureType2D = 2
   integer(c_int), parameter :: cudaTextureType2DLayered = 242
   integer(c_int), parameter :: cudaTextureType3D = 3
   integer(c_int), parameter :: cudaTextureTypeCubemap = 12
   integer(c_int), parameter :: cudaTextureTypeCubemapLayered = 252

   ! ======================================================================
   !  Interoperable derived types
   ! ======================================================================
   ! cudaLaunchAttributeValue is a C union: no Fortran equivalent, so it is
   ! declared as an opaque 64-byte buffer (size measured by the C probe).
   type, bind(C) :: cudaLaunchAttributeValue
      integer(c_int64_t) :: raw(8)
   end type cudaLaunchAttributeValue

   type, bind(C) :: char1
      integer(c_signed_char) :: x
   end type char1

   type, bind(C) :: uchar1
      integer(c_signed_char) :: x
   end type uchar1

   type, bind(C) :: char3
      integer(c_signed_char) :: x
      integer(c_signed_char) :: y
      integer(c_signed_char) :: z
   end type char3

   type, bind(C) :: uchar3
      integer(c_signed_char) :: x
      integer(c_signed_char) :: y
      integer(c_signed_char) :: z
   end type uchar3

   type, bind(C) :: short1
      integer(c_short) :: x
   end type short1

   type, bind(C) :: ushort1
      integer(c_short) :: x
   end type ushort1

   type, bind(C) :: short3
      integer(c_short) :: x
      integer(c_short) :: y
      integer(c_short) :: z
   end type short3

   type, bind(C) :: ushort3
      integer(c_short) :: x
      integer(c_short) :: y
      integer(c_short) :: z
   end type ushort3

   type, bind(C) :: int1
      integer(c_int) :: x
   end type int1

   type, bind(C) :: uint1
      integer(c_int) :: x
   end type uint1

   type, bind(C) :: int3
      integer(c_int) :: x
      integer(c_int) :: y
      integer(c_int) :: z
   end type int3

   type, bind(C) :: uint3
      integer(c_int) :: x
      integer(c_int) :: y
      integer(c_int) :: z
   end type uint3

   type, bind(C) :: long1
      integer(c_long) :: x
   end type long1

   type, bind(C) :: ulong1
      integer(c_long) :: x
   end type ulong1

   type, bind(C) :: long3
      integer(c_long) :: x
      integer(c_long) :: y
      integer(c_long) :: z
   end type long3

   type, bind(C) :: ulong3
      integer(c_long) :: x
      integer(c_long) :: y
      integer(c_long) :: z
   end type ulong3

   type, bind(C) :: float1
      real(c_float) :: x
   end type float1

   type, bind(C) :: float3
      real(c_float) :: x
      real(c_float) :: y
      real(c_float) :: z
   end type float3

   type, bind(C) :: longlong1
      integer(c_long_long) :: x
   end type longlong1

   type, bind(C) :: ulonglong1
      integer(c_long_long) :: x
   end type ulonglong1

   type, bind(C) :: longlong3
      integer(c_long_long) :: x
      integer(c_long_long) :: y
      integer(c_long_long) :: z
   end type longlong3

   type, bind(C) :: ulonglong3
      integer(c_long_long) :: x
      integer(c_long_long) :: y
      integer(c_long_long) :: z
   end type ulonglong3

   type, bind(C) :: double1
      real(c_double) :: x
   end type double1

   type, bind(C) :: double3
      real(c_double) :: x
      real(c_double) :: y
      real(c_double) :: z
   end type double3

   type, bind(C) :: dim3
      integer(c_int) :: x
      integer(c_int) :: y
      integer(c_int) :: z
   end type dim3

   type, bind(C) :: cudaChannelFormatDesc
      integer(c_int) :: x
      integer(c_int) :: y
      integer(c_int) :: z
      integer(c_int) :: w
      integer(c_int) :: f
   end type cudaChannelFormatDesc

   type, bind(C) :: cudaArraySparseProperties
      integer(c_int32_t) :: tileExtent(3)   ! C union / anonymous struct
      integer(c_int) :: miptailFirstLevel
      integer(c_long_long) :: miptailSize
      integer(c_int) :: flags
      integer(c_int) :: reserved(4)
   end type cudaArraySparseProperties

   type, bind(C) :: cudaArrayMemoryRequirements
      integer(c_size_t) :: size
      integer(c_size_t) :: alignment
      integer(c_int) :: reserved(4)
   end type cudaArrayMemoryRequirements

   type, bind(C) :: cudaPitchedPtr
      type(c_ptr) :: ptr
      integer(c_size_t) :: pitch
      integer(c_size_t) :: xsize
      integer(c_size_t) :: ysize
   end type cudaPitchedPtr

   type, bind(C) :: cudaExtent
      integer(c_size_t) :: width
      integer(c_size_t) :: height
      integer(c_size_t) :: depth
   end type cudaExtent

   type, bind(C) :: cudaPos
      integer(c_size_t) :: x
      integer(c_size_t) :: y
      integer(c_size_t) :: z
   end type cudaPos

   type, bind(C) :: cudaMemcpy3DParms
      type(c_ptr) :: srcArray
      type(cudaPos) :: srcPos
      type(cudaPitchedPtr) :: srcPtr
      type(c_ptr) :: dstArray
      type(cudaPos) :: dstPos
      type(cudaPitchedPtr) :: dstPtr
      type(cudaExtent) :: extent
      integer(c_int) :: kind
   end type cudaMemcpy3DParms

   type, bind(C) :: cudaMemcpyNodeParams
      integer(c_int) :: flags
      integer(c_int) :: reserved(3)
      type(cudaMemcpy3DParms) :: copyParams
   end type cudaMemcpyNodeParams

   type, bind(C) :: cudaMemcpy3DPeerParms
      type(c_ptr) :: srcArray
      type(cudaPos) :: srcPos
      type(cudaPitchedPtr) :: srcPtr
      integer(c_int) :: srcDevice
      type(c_ptr) :: dstArray
      type(cudaPos) :: dstPos
      type(cudaPitchedPtr) :: dstPtr
      integer(c_int) :: dstDevice
      type(cudaExtent) :: extent
   end type cudaMemcpy3DPeerParms

   type, bind(C) :: cudaMemsetParams
      type(c_ptr) :: dst
      integer(c_size_t) :: pitch
      integer(c_int) :: value
      integer(c_int) :: elementSize
      integer(c_size_t) :: width
      integer(c_size_t) :: height
   end type cudaMemsetParams

   type, bind(C) :: cudaMemsetParamsV2
      type(c_ptr) :: dst
      integer(c_size_t) :: pitch
      integer(c_int) :: value
      integer(c_int) :: elementSize
      integer(c_size_t) :: width
      integer(c_size_t) :: height
   end type cudaMemsetParamsV2

   type, bind(C) :: cudaAccessPolicyWindow
      type(c_ptr) :: base_ptr
      integer(c_size_t) :: num_bytes
      real(c_float) :: hitRatio
      integer(c_int) :: hitProp
      integer(c_int) :: missProp
   end type cudaAccessPolicyWindow

   type, bind(C) :: cudaHostNodeParams
      type(c_ptr) :: fn
      type(c_ptr) :: userData
   end type cudaHostNodeParams

   type, bind(C) :: cudaHostNodeParamsV2
      type(c_ptr) :: fn
      type(c_ptr) :: userData
   end type cudaHostNodeParamsV2

   type, bind(C) :: cudaResourceDesc
      integer(c_int) :: resType
      integer(c_int64_t) :: res(7)   ! C union / anonymous struct
   end type cudaResourceDesc

   type, bind(C) :: cudaResourceViewDesc
      integer(c_int) :: format
      integer(c_size_t) :: width
      integer(c_size_t) :: height
      integer(c_size_t) :: depth
      integer(c_int) :: firstMipmapLevel
      integer(c_int) :: lastMipmapLevel
      integer(c_int) :: firstLayer
      integer(c_int) :: lastLayer
   end type cudaResourceViewDesc

   type, bind(C) :: cudaPointerAttributes
      integer(c_int) :: type
      integer(c_int) :: device
      type(c_ptr) :: devicePointer
      type(c_ptr) :: hostPointer
   end type cudaPointerAttributes

   type, bind(C) :: cudaFuncAttributes
      integer(c_size_t) :: sharedSizeBytes
      integer(c_size_t) :: constSizeBytes
      integer(c_size_t) :: localSizeBytes
      integer(c_int) :: maxThreadsPerBlock
      integer(c_int) :: numRegs
      integer(c_int) :: ptxVersion
      integer(c_int) :: binaryVersion
      integer(c_int) :: cacheModeCA
      integer(c_int) :: maxDynamicSharedSizeBytes
      integer(c_int) :: preferredShmemCarveout
      integer(c_int) :: clusterDimMustBeSet
      integer(c_int) :: requiredClusterWidth
      integer(c_int) :: requiredClusterHeight
      integer(c_int) :: requiredClusterDepth
      integer(c_int) :: clusterSchedulingPolicyPreference
      integer(c_int) :: nonPortableClusterSizeAllowed
      integer(c_int) :: reserved(16)
   end type cudaFuncAttributes

   type, bind(C) :: cudaMemLocation
      integer(c_int) :: type
      integer(c_int) :: id
   end type cudaMemLocation

   type, bind(C) :: cudaMemAccessDesc
      type(cudaMemLocation) :: location
      integer(c_int) :: flags
   end type cudaMemAccessDesc

   type, bind(C) :: cudaMemPoolProps
      integer(c_int) :: allocType
      integer(c_int) :: handleTypes
      type(cudaMemLocation) :: location
      type(c_ptr) :: win32SecurityAttributes
      integer(c_size_t) :: maxSize
      integer(c_short) :: usage
      integer(c_signed_char) :: reserved(54)
   end type cudaMemPoolProps

   type, bind(C) :: cudaMemPoolPtrExportData
      integer(c_signed_char) :: reserved(64)
   end type cudaMemPoolPtrExportData

   type, bind(C) :: cudaMemAllocNodeParams
      type(cudaMemPoolProps) :: poolProps
      type(c_ptr) :: accessDescs
      integer(c_size_t) :: accessDescCount
      integer(c_size_t) :: bytesize
      type(c_ptr) :: dptr
   end type cudaMemAllocNodeParams

   type, bind(C) :: cudaMemAllocNodeParamsV2
      type(cudaMemPoolProps) :: poolProps
      type(c_ptr) :: accessDescs
      integer(c_size_t) :: accessDescCount
      integer(c_size_t) :: bytesize
      type(c_ptr) :: dptr
   end type cudaMemAllocNodeParamsV2

   type, bind(C) :: cudaMemFreeNodeParams
      type(c_ptr) :: dptr
   end type cudaMemFreeNodeParams

   type, bind(C) :: cudaMemcpyAttributes
      integer(c_int) :: srcAccessOrder
      type(cudaMemLocation) :: srcLocHint
      type(cudaMemLocation) :: dstLocHint
      integer(c_int) :: flags
   end type cudaMemcpyAttributes

   type, bind(C) :: cudaOffset3D
      integer(c_size_t) :: x
      integer(c_size_t) :: y
      integer(c_size_t) :: z
   end type cudaOffset3D

   type, bind(C) :: cudaMemcpy3DOperand
      integer(c_int) :: type
      integer(c_int64_t) :: op(4)   ! C union / anonymous struct
   end type cudaMemcpy3DOperand

   type, bind(C) :: cudaMemcpy3DBatchOp
      type(cudaMemcpy3DOperand) :: src
      type(cudaMemcpy3DOperand) :: dst
      type(cudaExtent) :: extent
      integer(c_int) :: srcAccessOrder
      integer(c_int) :: flags
   end type cudaMemcpy3DBatchOp

   type, bind(C) :: CUuuid_st
      character(kind=c_char) :: bytes(16)
   end type CUuuid_st

   type, bind(C) :: cudaDeviceProp
      character(kind=c_char) :: name(256)
      type(CUuuid_st) :: uuid
      character(kind=c_char) :: luid(8)
      integer(c_int) :: luidDeviceNodeMask
      integer(c_size_t) :: totalGlobalMem
      integer(c_size_t) :: sharedMemPerBlock
      integer(c_int) :: regsPerBlock
      integer(c_int) :: warpSize
      integer(c_size_t) :: memPitch
      integer(c_int) :: maxThreadsPerBlock
      integer(c_int) :: maxThreadsDim(3)
      integer(c_int) :: maxGridSize(3)
      integer(c_int) :: clockRate
      integer(c_size_t) :: totalConstMem
      integer(c_int) :: major
      integer(c_int) :: minor
      integer(c_size_t) :: textureAlignment
      integer(c_size_t) :: texturePitchAlignment
      integer(c_int) :: deviceOverlap
      integer(c_int) :: multiProcessorCount
      integer(c_int) :: kernelExecTimeoutEnabled
      integer(c_int) :: integrated
      integer(c_int) :: canMapHostMemory
      integer(c_int) :: computeMode
      integer(c_int) :: maxTexture1D
      integer(c_int) :: maxTexture1DMipmap
      integer(c_int) :: maxTexture1DLinear
      integer(c_int) :: maxTexture2D(2)
      integer(c_int) :: maxTexture2DMipmap(2)
      integer(c_int) :: maxTexture2DLinear(3)
      integer(c_int) :: maxTexture2DGather(2)
      integer(c_int) :: maxTexture3D(3)
      integer(c_int) :: maxTexture3DAlt(3)
      integer(c_int) :: maxTextureCubemap
      integer(c_int) :: maxTexture1DLayered(2)
      integer(c_int) :: maxTexture2DLayered(3)
      integer(c_int) :: maxTextureCubemapLayered(2)
      integer(c_int) :: maxSurface1D
      integer(c_int) :: maxSurface2D(2)
      integer(c_int) :: maxSurface3D(3)
      integer(c_int) :: maxSurface1DLayered(2)
      integer(c_int) :: maxSurface2DLayered(3)
      integer(c_int) :: maxSurfaceCubemap
      integer(c_int) :: maxSurfaceCubemapLayered(2)
      integer(c_size_t) :: surfaceAlignment
      integer(c_int) :: concurrentKernels
      integer(c_int) :: ECCEnabled
      integer(c_int) :: pciBusID
      integer(c_int) :: pciDeviceID
      integer(c_int) :: pciDomainID
      integer(c_int) :: tccDriver
      integer(c_int) :: asyncEngineCount
      integer(c_int) :: unifiedAddressing
      integer(c_int) :: memoryClockRate
      integer(c_int) :: memoryBusWidth
      integer(c_int) :: l2CacheSize
      integer(c_int) :: persistingL2CacheMaxSize
      integer(c_int) :: maxThreadsPerMultiProcessor
      integer(c_int) :: streamPrioritiesSupported
      integer(c_int) :: globalL1CacheSupported
      integer(c_int) :: localL1CacheSupported
      integer(c_size_t) :: sharedMemPerMultiprocessor
      integer(c_int) :: regsPerMultiprocessor
      integer(c_int) :: managedMemory
      integer(c_int) :: isMultiGpuBoard
      integer(c_int) :: multiGpuBoardGroupID
      integer(c_int) :: hostNativeAtomicSupported
      integer(c_int) :: singleToDoublePrecisionPerfRatio
      integer(c_int) :: pageableMemoryAccess
      integer(c_int) :: concurrentManagedAccess
      integer(c_int) :: computePreemptionSupported
      integer(c_int) :: canUseHostPointerForRegisteredMem
      integer(c_int) :: cooperativeLaunch
      integer(c_int) :: cooperativeMultiDeviceLaunch
      integer(c_size_t) :: sharedMemPerBlockOptin
      integer(c_int) :: pageableMemoryAccessUsesHostPageTables
      integer(c_int) :: directManagedMemAccessFromHost
      integer(c_int) :: maxBlocksPerMultiProcessor
      integer(c_int) :: accessPolicyMaxWindowSize
      integer(c_size_t) :: reservedSharedMemPerBlock
      integer(c_int) :: hostRegisterSupported
      integer(c_int) :: sparseCudaArraySupported
      integer(c_int) :: hostRegisterReadOnlySupported
      integer(c_int) :: timelineSemaphoreInteropSupported
      integer(c_int) :: memoryPoolsSupported
      integer(c_int) :: gpuDirectRDMASupported
      integer(c_int) :: gpuDirectRDMAFlushWritesOptions
      integer(c_int) :: gpuDirectRDMAWritesOrdering
      integer(c_int) :: memoryPoolSupportedHandleTypes
      integer(c_int) :: deferredMappingCudaArraySupported
      integer(c_int) :: ipcEventSupported
      integer(c_int) :: clusterLaunch
      integer(c_int) :: unifiedFunctionPointers
      integer(c_int) :: reserved(63)
   end type cudaDeviceProp

   type, bind(C) :: cudaIpcEventHandle_st
      character(kind=c_char) :: reserved(64)
   end type cudaIpcEventHandle_st

   type, bind(C) :: cudaIpcMemHandle_st
      character(kind=c_char) :: reserved(64)
   end type cudaIpcMemHandle_st

   type, bind(C) :: cudaMemFabricHandle_st
      character(kind=c_char) :: reserved(64)
   end type cudaMemFabricHandle_st

   type, bind(C) :: cudaExternalMemoryHandleDesc
      integer(c_int) :: type
      integer(c_int64_t) :: handle(2)   ! C union / anonymous struct
      integer(c_long_long) :: size
      integer(c_int) :: flags
   end type cudaExternalMemoryHandleDesc

   type, bind(C) :: cudaExternalMemoryBufferDesc
      integer(c_long_long) :: offset
      integer(c_long_long) :: size
      integer(c_int) :: flags
   end type cudaExternalMemoryBufferDesc

   type, bind(C) :: cudaExternalMemoryMipmappedArrayDesc
      integer(c_long_long) :: offset
      type(cudaChannelFormatDesc) :: formatDesc
      type(cudaExtent) :: extent
      integer(c_int) :: flags
      integer(c_int) :: numLevels
   end type cudaExternalMemoryMipmappedArrayDesc

   type, bind(C) :: cudaExternalSemaphoreHandleDesc
      integer(c_int) :: type
      integer(c_int64_t) :: handle(2)   ! C union / anonymous struct
      integer(c_int) :: flags
   end type cudaExternalSemaphoreHandleDesc

   type, bind(C) :: cudaExternalSemaphoreSignalParams_v1
      integer(c_int64_t) :: params(3)   ! C union / anonymous struct
      integer(c_int) :: flags
   end type cudaExternalSemaphoreSignalParams_v1

   type, bind(C) :: cudaExternalSemaphoreWaitParams_v1
      integer(c_int64_t) :: params(4)   ! C union / anonymous struct
      integer(c_int) :: flags
   end type cudaExternalSemaphoreWaitParams_v1

   type, bind(C) :: cudaExternalSemaphoreSignalParams
      integer(c_int64_t) :: params(9)   ! C union / anonymous struct
      integer(c_int) :: flags
      integer(c_int) :: reserved(16)
   end type cudaExternalSemaphoreSignalParams

   type, bind(C) :: cudaExternalSemaphoreWaitParams
      integer(c_int64_t) :: params(9)   ! C union / anonymous struct
      integer(c_int) :: flags
      integer(c_int) :: reserved(16)
   end type cudaExternalSemaphoreWaitParams

   type, bind(C) :: cudaLaunchParams
      type(c_ptr) :: func
      type(dim3) :: gridDim
      type(dim3) :: blockDim
      type(c_ptr) :: args
      integer(c_size_t) :: sharedMem
      type(c_ptr) :: stream
   end type cudaLaunchParams

   type, bind(C) :: cudaKernelNodeParams
      type(c_ptr) :: func
      type(dim3) :: gridDim
      type(dim3) :: blockDim
      integer(c_int) :: sharedMemBytes
      type(c_ptr) :: kernelParams
      type(c_ptr) :: extra
   end type cudaKernelNodeParams

   type, bind(C) :: cudaKernelNodeParamsV2
      type(c_ptr) :: func
      type(dim3) :: gridDim
      type(dim3) :: blockDim
      integer(c_int) :: sharedMemBytes
      type(c_ptr) :: kernelParams
      type(c_ptr) :: extra
   end type cudaKernelNodeParamsV2

   type, bind(C) :: cudaExternalSemaphoreSignalNodeParams
      type(c_ptr) :: extSemArray
      type(c_ptr) :: paramsArray
      integer(c_int) :: numExtSems
   end type cudaExternalSemaphoreSignalNodeParams

   type, bind(C) :: cudaExternalSemaphoreSignalNodeParamsV2
      type(c_ptr) :: extSemArray
      type(c_ptr) :: paramsArray
      integer(c_int) :: numExtSems
   end type cudaExternalSemaphoreSignalNodeParamsV2

   type, bind(C) :: cudaExternalSemaphoreWaitNodeParams
      type(c_ptr) :: extSemArray
      type(c_ptr) :: paramsArray
      integer(c_int) :: numExtSems
   end type cudaExternalSemaphoreWaitNodeParams

   type, bind(C) :: cudaExternalSemaphoreWaitNodeParamsV2
      type(c_ptr) :: extSemArray
      type(c_ptr) :: paramsArray
      integer(c_int) :: numExtSems
   end type cudaExternalSemaphoreWaitNodeParamsV2

   type, bind(C) :: cudaConditionalNodeParams
      integer(c_long_long) :: handle
      integer(c_int) :: type
      integer(c_int) :: size
      type(c_ptr) :: phGraph_out
   end type cudaConditionalNodeParams

   type, bind(C) :: cudaChildGraphNodeParams
      type(c_ptr) :: graph
      integer(c_int) :: ownership
   end type cudaChildGraphNodeParams

   type, bind(C) :: cudaEventRecordNodeParams
      type(c_ptr) :: event
   end type cudaEventRecordNodeParams

   type, bind(C) :: cudaEventWaitNodeParams
      type(c_ptr) :: event
   end type cudaEventWaitNodeParams

   type, bind(C) :: cudaGraphNodeParams
      integer(c_int) :: type
      integer(c_int) :: reserved0(3)
      integer(c_int8_t) :: anon2(232)   ! anonymous C union (+ padding)
      integer(c_long_long) :: reserved2
   end type cudaGraphNodeParams

   type, bind(C) :: cudaGraphEdgeData_st
      integer(c_signed_char) :: from_port
      integer(c_signed_char) :: to_port
      integer(c_signed_char) :: type
      integer(c_signed_char) :: reserved(5)
   end type cudaGraphEdgeData_st

   type, bind(C) :: cudaGraphInstantiateParams_st
      integer(c_long_long) :: flags
      type(c_ptr) :: uploadStream
      type(c_ptr) :: errNode_out
      integer(c_int) :: result_out
   end type cudaGraphInstantiateParams_st

   type, bind(C) :: cudaGraphExecUpdateResultInfo_st
      integer(c_int) :: result
      type(c_ptr) :: errorNode
      type(c_ptr) :: errorFromNode
   end type cudaGraphExecUpdateResultInfo_st

   type, bind(C) :: cudaGraphKernelNodeUpdate
      type(c_ptr) :: node
      integer(c_int) :: field
      integer(c_int64_t) :: updateData(3)   ! C union / anonymous struct
   end type cudaGraphKernelNodeUpdate

   type, bind(C) :: cudaLaunchMemSyncDomainMap_st
      integer(c_signed_char) :: default_
      integer(c_signed_char) :: remote
   end type cudaLaunchMemSyncDomainMap_st

   type, bind(C) :: cudaLaunchAttribute_st
      integer(c_int) :: id
      integer(c_int8_t) :: pad(4)   ! computed-size C array (+ padding)
      type(cudaLaunchAttributeValue) :: val
   end type cudaLaunchAttribute_st

   type, bind(C) :: cudaLaunchConfig_st
      type(dim3) :: gridDim
      type(dim3) :: blockDim
      integer(c_size_t) :: dynamicSmemBytes
      type(c_ptr) :: stream
      type(c_ptr) :: attrs
      integer(c_int) :: numAttrs
   end type cudaLaunchConfig_st

   type, bind(C) :: cudaAsyncNotificationInfo
      integer(c_int) :: type
      integer(c_int64_t) :: info(1)   ! C union / anonymous struct
   end type cudaAsyncNotificationInfo

   type, bind(C) :: cudaTextureDesc
      integer(c_int) :: addressMode(3)
      integer(c_int) :: filterMode
      integer(c_int) :: readMode
      integer(c_int) :: sRGB
      real(c_float) :: borderColor(4)
      integer(c_int) :: normalizedCoords
      integer(c_int) :: maxAnisotropy
      integer(c_int) :: mipmapFilterMode
      real(c_float) :: mipmapLevelBias
      real(c_float) :: minMipmapLevelClamp
      real(c_float) :: maxMipmapLevelClamp
      integer(c_int) :: disableTrilinearOptimization
      integer(c_int) :: seamlessCubemap
   end type cudaTextureDesc

   type, bind(C) :: cudalibraryHostUniversalFunctionAndDataTable_t
      type(c_ptr) :: functionTable
      integer(c_size_t) :: functionWindowSize
      type(c_ptr) :: dataTable
      integer(c_size_t) :: dataWindowSize
   end type cudalibraryHostUniversalFunctionAndDataTable_t

   ! ======================================================================
   !  C entry points
   ! ======================================================================
   interface

      integer(c_int) function cudaArrayGetInfo(desc, extent, flags, array) &
         bind(C, name="cudaArrayGetInfo")
         import
         type(cudaChannelFormatDesc), intent(inout) :: desc
         type(cudaExtent), intent(inout) :: extent
         integer(c_int), intent(inout) :: flags
         type(c_ptr), value :: array
      end function cudaArrayGetInfo

      integer(c_int) function cudaArrayGetMemoryRequirements(memoryRequirements, array, device) &
         bind(C, name="cudaArrayGetMemoryRequirements")
         import
         type(cudaArrayMemoryRequirements), intent(inout) :: memoryRequirements
         type(c_ptr), value :: array
         integer(c_int), value :: device
      end function cudaArrayGetMemoryRequirements

      integer(c_int) function cudaArrayGetPlane(pPlaneArray, hArray, planeIdx) &
         bind(C, name="cudaArrayGetPlane")
         import
         type(c_ptr), intent(out) :: pPlaneArray
         type(c_ptr), value :: hArray
         integer(c_int), value :: planeIdx
      end function cudaArrayGetPlane

      integer(c_int) function cudaArrayGetSparseProperties(sparseProperties, array) &
         bind(C, name="cudaArrayGetSparseProperties")
         import
         type(cudaArraySparseProperties), intent(inout) :: sparseProperties
         type(c_ptr), value :: array
      end function cudaArrayGetSparseProperties

      integer(c_int) function cudaChooseDevice(device, prop) &
         bind(C, name="cudaChooseDevice")
         import
         integer(c_int), intent(inout) :: device
         type(cudaDeviceProp), intent(in) :: prop
      end function cudaChooseDevice

      integer(c_int) function cudaCreateSurfaceObject(pSurfObject, pResDesc) &
         bind(C, name="cudaCreateSurfaceObject")
         import
         integer(c_long_long), intent(inout) :: pSurfObject
         type(cudaResourceDesc), intent(in) :: pResDesc
      end function cudaCreateSurfaceObject

      integer(c_int) function cudaCreateTextureObject(pTexObject, pResDesc, pTexDesc, pResViewDesc) &
         bind(C, name="cudaCreateTextureObject")
         import
         integer(c_long_long), intent(inout) :: pTexObject
         type(cudaResourceDesc), intent(in) :: pResDesc
         type(cudaTextureDesc), intent(in) :: pTexDesc
         type(cudaResourceViewDesc), intent(in) :: pResViewDesc
      end function cudaCreateTextureObject

      integer(c_int) function cudaCtxResetPersistingL2Cache() &
         bind(C, name="cudaCtxResetPersistingL2Cache")
         import
      end function cudaCtxResetPersistingL2Cache

      integer(c_int) function cudaDestroyExternalMemory(extMem) &
         bind(C, name="cudaDestroyExternalMemory")
         import
         type(c_ptr), value :: extMem
      end function cudaDestroyExternalMemory

      integer(c_int) function cudaDestroyExternalSemaphore(extSem) &
         bind(C, name="cudaDestroyExternalSemaphore")
         import
         type(c_ptr), value :: extSem
      end function cudaDestroyExternalSemaphore

      integer(c_int) function cudaDestroySurfaceObject(surfObject) &
         bind(C, name="cudaDestroySurfaceObject")
         import
         integer(c_long_long), value :: surfObject
      end function cudaDestroySurfaceObject

      integer(c_int) function cudaDestroyTextureObject(texObject) &
         bind(C, name="cudaDestroyTextureObject")
         import
         integer(c_long_long), value :: texObject
      end function cudaDestroyTextureObject

      integer(c_int) function cudaDeviceCanAccessPeer(canAccessPeer, device, peerDevice) &
         bind(C, name="cudaDeviceCanAccessPeer")
         import
         integer(c_int), intent(inout) :: canAccessPeer
         integer(c_int), value :: device
         integer(c_int), value :: peerDevice
      end function cudaDeviceCanAccessPeer

      integer(c_int) function cudaDeviceDisablePeerAccess(peerDevice) &
         bind(C, name="cudaDeviceDisablePeerAccess")
         import
         integer(c_int), value :: peerDevice
      end function cudaDeviceDisablePeerAccess

      integer(c_int) function cudaDeviceEnablePeerAccess(peerDevice, flags) &
         bind(C, name="cudaDeviceEnablePeerAccess")
         import
         integer(c_int), value :: peerDevice
         integer(c_int), value :: flags
      end function cudaDeviceEnablePeerAccess

      integer(c_int) function cudaDeviceFlushGPUDirectRDMAWrites(target, scope) &
         bind(C, name="cudaDeviceFlushGPUDirectRDMAWrites")
         import
         integer(c_int), value :: target
         integer(c_int), value :: scope
      end function cudaDeviceFlushGPUDirectRDMAWrites

      integer(c_int) function cudaDeviceGetAttribute(value, attr, device) &
         bind(C, name="cudaDeviceGetAttribute")
         import
         integer(c_int), intent(inout) :: value
         integer(c_int), value :: attr
         integer(c_int), value :: device
      end function cudaDeviceGetAttribute

      integer(c_int) function cudaDeviceGetByPCIBusId(device, pciBusId) &
         bind(C, name="cudaDeviceGetByPCIBusId")
         import
         integer(c_int), intent(inout) :: device
         character(kind=c_char), dimension(*), intent(in) :: pciBusId
      end function cudaDeviceGetByPCIBusId

      integer(c_int) function cudaDeviceGetCacheConfig(pCacheConfig) &
         bind(C, name="cudaDeviceGetCacheConfig")
         import
         integer(c_int), intent(out) :: pCacheConfig
      end function cudaDeviceGetCacheConfig

      integer(c_int) function cudaDeviceGetDefaultMemPool(memPool, device) &
         bind(C, name="cudaDeviceGetDefaultMemPool")
         import
         type(c_ptr), intent(out) :: memPool
         integer(c_int), value :: device
      end function cudaDeviceGetDefaultMemPool

      integer(c_int) function cudaDeviceGetGraphMemAttribute(device, attr, value) &
         bind(C, name="cudaDeviceGetGraphMemAttribute")
         import
         integer(c_int), value :: device
         integer(c_int), value :: attr
         type(c_ptr), value :: value
      end function cudaDeviceGetGraphMemAttribute

      integer(c_int) function cudaDeviceGetLimit(pValue, limit) &
         bind(C, name="cudaDeviceGetLimit")
         import
         integer(c_size_t), intent(inout) :: pValue
         integer(c_int), value :: limit
      end function cudaDeviceGetLimit

      integer(c_int) function cudaDeviceGetMemPool(memPool, device) &
         bind(C, name="cudaDeviceGetMemPool")
         import
         type(c_ptr), intent(out) :: memPool
         integer(c_int), value :: device
      end function cudaDeviceGetMemPool

      integer(c_int) function cudaDeviceGetNvSciSyncAttributes(nvSciSyncAttrList, device, flags) &
         bind(C, name="cudaDeviceGetNvSciSyncAttributes")
         import
         type(c_ptr), value :: nvSciSyncAttrList
         integer(c_int), value :: device
         integer(c_int), value :: flags
      end function cudaDeviceGetNvSciSyncAttributes

      integer(c_int) function cudaDeviceGetP2PAttribute(value, attr, srcDevice, dstDevice) &
         bind(C, name="cudaDeviceGetP2PAttribute")
         import
         integer(c_int), intent(inout) :: value
         integer(c_int), value :: attr
         integer(c_int), value :: srcDevice
         integer(c_int), value :: dstDevice
      end function cudaDeviceGetP2PAttribute

      integer(c_int) function cudaDeviceGetPCIBusId(pciBusId, len, device) &
         bind(C, name="cudaDeviceGetPCIBusId")
         import
         character(kind=c_char), dimension(*), intent(in) :: pciBusId
         integer(c_int), value :: len
         integer(c_int), value :: device
      end function cudaDeviceGetPCIBusId

      integer(c_int) function cudaDeviceGetSharedMemConfig(pConfig) &
         bind(C, name="cudaDeviceGetSharedMemConfig")
         import
         integer(c_int), intent(out) :: pConfig
      end function cudaDeviceGetSharedMemConfig

      integer(c_int) function cudaDeviceGetStreamPriorityRange(leastPriority, greatestPriority) &
         bind(C, name="cudaDeviceGetStreamPriorityRange")
         import
         integer(c_int), intent(inout) :: leastPriority
         integer(c_int), intent(inout) :: greatestPriority
      end function cudaDeviceGetStreamPriorityRange

      integer(c_int) function cudaDeviceGetTexture1DLinearMaxWidth(maxWidthInElements, fmtDesc, device) &
         bind(C, name="cudaDeviceGetTexture1DLinearMaxWidth")
         import
         integer(c_size_t), intent(inout) :: maxWidthInElements
         type(cudaChannelFormatDesc), intent(in) :: fmtDesc
         integer(c_int), value :: device
      end function cudaDeviceGetTexture1DLinearMaxWidth

      integer(c_int) function cudaDeviceGraphMemTrim(device) &
         bind(C, name="cudaDeviceGraphMemTrim")
         import
         integer(c_int), value :: device
      end function cudaDeviceGraphMemTrim

      integer(c_int) function cudaDeviceRegisterAsyncNotification(device, callbackFunc, userData, callback) &
         bind(C, name="cudaDeviceRegisterAsyncNotification")
         import
         integer(c_int), value :: device
         type(c_funptr), value :: callbackFunc
         type(c_ptr), value :: userData
         type(c_ptr), intent(out) :: callback
      end function cudaDeviceRegisterAsyncNotification

      integer(c_int) function cudaDeviceReset() &
         bind(C, name="cudaDeviceReset")
         import
      end function cudaDeviceReset

      integer(c_int) function cudaDeviceSetCacheConfig(cacheConfig) &
         bind(C, name="cudaDeviceSetCacheConfig")
         import
         integer(c_int), value :: cacheConfig
      end function cudaDeviceSetCacheConfig

      integer(c_int) function cudaDeviceSetGraphMemAttribute(device, attr, value) &
         bind(C, name="cudaDeviceSetGraphMemAttribute")
         import
         integer(c_int), value :: device
         integer(c_int), value :: attr
         type(c_ptr), value :: value
      end function cudaDeviceSetGraphMemAttribute

      integer(c_int) function cudaDeviceSetLimit(limit, value) &
         bind(C, name="cudaDeviceSetLimit")
         import
         integer(c_int), value :: limit
         integer(c_size_t), value :: value
      end function cudaDeviceSetLimit

      integer(c_int) function cudaDeviceSetMemPool(device, memPool) &
         bind(C, name="cudaDeviceSetMemPool")
         import
         integer(c_int), value :: device
         type(c_ptr), value :: memPool
      end function cudaDeviceSetMemPool

      integer(c_int) function cudaDeviceSetSharedMemConfig(config) &
         bind(C, name="cudaDeviceSetSharedMemConfig")
         import
         integer(c_int), value :: config
      end function cudaDeviceSetSharedMemConfig

      integer(c_int) function cudaDeviceSynchronize() &
         bind(C, name="cudaDeviceSynchronize")
         import
      end function cudaDeviceSynchronize

      integer(c_int) function cudaDeviceUnregisterAsyncNotification(device, callback) &
         bind(C, name="cudaDeviceUnregisterAsyncNotification")
         import
         integer(c_int), value :: device
         type(c_ptr), value :: callback
      end function cudaDeviceUnregisterAsyncNotification

      integer(c_int) function cudaDriverGetVersion(driverVersion) &
         bind(C, name="cudaDriverGetVersion")
         import
         integer(c_int), intent(inout) :: driverVersion
      end function cudaDriverGetVersion

      integer(c_int) function cudaEventCreate(event) &
         bind(C, name="cudaEventCreate")
         import
         type(c_ptr), intent(out) :: event
      end function cudaEventCreate

      integer(c_int) function cudaEventCreateWithFlags(event, flags) &
         bind(C, name="cudaEventCreateWithFlags")
         import
         type(c_ptr), intent(out) :: event
         integer(c_int), value :: flags
      end function cudaEventCreateWithFlags

      integer(c_int) function cudaEventDestroy(event) &
         bind(C, name="cudaEventDestroy")
         import
         type(c_ptr), value :: event
      end function cudaEventDestroy

      integer(c_int) function cudaEventElapsedTime(ms, start, end) &
         bind(C, name="cudaEventElapsedTime")
         import
         real(c_float), intent(inout) :: ms
         type(c_ptr), value :: start
         type(c_ptr), value :: end
      end function cudaEventElapsedTime

      integer(c_int) function cudaEventElapsedTime_v2(ms, start, end) &
         bind(C, name="cudaEventElapsedTime_v2")
         import
         real(c_float), intent(inout) :: ms
         type(c_ptr), value :: start
         type(c_ptr), value :: end
      end function cudaEventElapsedTime_v2

      integer(c_int) function cudaEventQuery(event) &
         bind(C, name="cudaEventQuery")
         import
         type(c_ptr), value :: event
      end function cudaEventQuery

      integer(c_int) function cudaEventRecord(event, stream) &
         bind(C, name="cudaEventRecord")
         import
         type(c_ptr), value :: event
         type(c_ptr), value :: stream
      end function cudaEventRecord

      integer(c_int) function cudaEventRecordWithFlags(event, stream, flags) &
         bind(C, name="cudaEventRecordWithFlags")
         import
         type(c_ptr), value :: event
         type(c_ptr), value :: stream
         integer(c_int), value :: flags
      end function cudaEventRecordWithFlags

      integer(c_int) function cudaEventSynchronize(event) &
         bind(C, name="cudaEventSynchronize")
         import
         type(c_ptr), value :: event
      end function cudaEventSynchronize

      integer(c_int) function cudaExternalMemoryGetMappedBuffer(devPtr, extMem, bufferDesc) &
         bind(C, name="cudaExternalMemoryGetMappedBuffer")
         import
         type(c_ptr), intent(out) :: devPtr
         type(c_ptr), value :: extMem
         type(cudaExternalMemoryBufferDesc), intent(in) :: bufferDesc
      end function cudaExternalMemoryGetMappedBuffer

      integer(c_int) function cudaExternalMemoryGetMappedMipmappedArray(mipmap, extMem, mipmapDesc) &
         bind(C, name="cudaExternalMemoryGetMappedMipmappedArray")
         import
         type(c_ptr), intent(out) :: mipmap
         type(c_ptr), value :: extMem
         type(cudaExternalMemoryMipmappedArrayDesc), intent(in) :: mipmapDesc
      end function cudaExternalMemoryGetMappedMipmappedArray

      integer(c_int) function cudaFree(devPtr) &
         bind(C, name="cudaFree")
         import
         type(c_ptr), value :: devPtr
      end function cudaFree

      integer(c_int) function cudaFreeArray(array) &
         bind(C, name="cudaFreeArray")
         import
         type(c_ptr), value :: array
      end function cudaFreeArray

      integer(c_int) function cudaFreeAsync(devPtr, hStream) &
         bind(C, name="cudaFreeAsync")
         import
         type(c_ptr), value :: devPtr
         type(c_ptr), value :: hStream
      end function cudaFreeAsync

      integer(c_int) function cudaFreeHost(ptr) &
         bind(C, name="cudaFreeHost")
         import
         type(c_ptr), value :: ptr
      end function cudaFreeHost

      integer(c_int) function cudaFreeMipmappedArray(mipmappedArray) &
         bind(C, name="cudaFreeMipmappedArray")
         import
         type(c_ptr), value :: mipmappedArray
      end function cudaFreeMipmappedArray

      integer(c_int) function cudaFuncGetAttributes(attr, func) &
         bind(C, name="cudaFuncGetAttributes")
         import
         type(cudaFuncAttributes), intent(inout) :: attr
         type(c_ptr), value :: func
      end function cudaFuncGetAttributes

      integer(c_int) function cudaFuncGetName(name, func) &
         bind(C, name="cudaFuncGetName")
         import
         type(c_ptr), intent(out) :: name
         type(c_ptr), value :: func
      end function cudaFuncGetName

      integer(c_int) function cudaFuncGetParamInfo(func, paramIndex, paramOffset, paramSize) &
         bind(C, name="cudaFuncGetParamInfo")
         import
         type(c_ptr), value :: func
         integer(c_size_t), value :: paramIndex
         integer(c_size_t), intent(inout) :: paramOffset
         integer(c_size_t), intent(inout) :: paramSize
      end function cudaFuncGetParamInfo

      integer(c_int) function cudaFuncSetAttribute(func, attr, value) &
         bind(C, name="cudaFuncSetAttribute")
         import
         type(c_ptr), value :: func
         integer(c_int), value :: attr
         integer(c_int), value :: value
      end function cudaFuncSetAttribute

      integer(c_int) function cudaFuncSetCacheConfig(func, cacheConfig) &
         bind(C, name="cudaFuncSetCacheConfig")
         import
         type(c_ptr), value :: func
         integer(c_int), value :: cacheConfig
      end function cudaFuncSetCacheConfig

      integer(c_int) function cudaFuncSetSharedMemConfig(func, config) &
         bind(C, name="cudaFuncSetSharedMemConfig")
         import
         type(c_ptr), value :: func
         integer(c_int), value :: config
      end function cudaFuncSetSharedMemConfig

      integer(c_int) function cudaGetChannelDesc(desc, array) &
         bind(C, name="cudaGetChannelDesc")
         import
         type(cudaChannelFormatDesc), intent(inout) :: desc
         type(c_ptr), value :: array
      end function cudaGetChannelDesc

      integer(c_int) function cudaGetDevice(device) &
         bind(C, name="cudaGetDevice")
         import
         integer(c_int), intent(inout) :: device
      end function cudaGetDevice

      integer(c_int) function cudaGetDeviceCount(count) &
         bind(C, name="cudaGetDeviceCount")
         import
         integer(c_int), intent(inout) :: count
      end function cudaGetDeviceCount

      integer(c_int) function cudaGetDeviceFlags(flags) &
         bind(C, name="cudaGetDeviceFlags")
         import
         integer(c_int), intent(inout) :: flags
      end function cudaGetDeviceFlags

      integer(c_int) function cudaGetDeviceProperties(prop, device) &
         bind(C, name="cudaGetDeviceProperties_v2")   ! header aliases cudaGetDeviceProperties -> cudaGetDeviceProperties_v2
         import
         type(cudaDeviceProp), intent(inout) :: prop
         integer(c_int), value :: device
      end function cudaGetDeviceProperties

      integer(c_int) function cudaGetDriverEntryPoint(symbol, funcPtr, flags, driverStatus) &
         bind(C, name="cudaGetDriverEntryPoint")
         import
         character(kind=c_char), dimension(*), intent(in) :: symbol
         type(c_ptr), intent(out) :: funcPtr
         integer(c_long_long), value :: flags
         integer(c_int), intent(out) :: driverStatus
      end function cudaGetDriverEntryPoint

      integer(c_int) function cudaGetDriverEntryPointByVersion(symbol, funcPtr, cudaVersion, flags, driverStatus) &
         bind(C, name="cudaGetDriverEntryPointByVersion")
         import
         character(kind=c_char), dimension(*), intent(in) :: symbol
         type(c_ptr), intent(out) :: funcPtr
         integer(c_int), value :: cudaVersion
         integer(c_long_long), value :: flags
         integer(c_int), intent(out) :: driverStatus
      end function cudaGetDriverEntryPointByVersion

      type(c_ptr) function cudaGetErrorName(error) &
         bind(C, name="cudaGetErrorName")
         import
         integer(c_int), value :: error
      end function cudaGetErrorName

      type(c_ptr) function cudaGetErrorString(error) &
         bind(C, name="cudaGetErrorString")
         import
         integer(c_int), value :: error
      end function cudaGetErrorString

      integer(c_int) function cudaGetExportTable(ppExportTable, pExportTableId) &
         bind(C, name="cudaGetExportTable")
         import
         type(c_ptr), intent(out) :: ppExportTable
         type(CUuuid_st), intent(in) :: pExportTableId
      end function cudaGetExportTable

      integer(c_int) function cudaGetFuncBySymbol(functionPtr, symbolPtr) &
         bind(C, name="cudaGetFuncBySymbol")
         import
         type(c_ptr), intent(out) :: functionPtr
         type(c_ptr), value :: symbolPtr
      end function cudaGetFuncBySymbol

      integer(c_int) function cudaGetKernel(kernelPtr, entryFuncAddr) &
         bind(C, name="cudaGetKernel")
         import
         type(c_ptr), intent(out) :: kernelPtr
         type(c_ptr), value :: entryFuncAddr
      end function cudaGetKernel

      integer(c_int) function cudaGetLastError() &
         bind(C, name="cudaGetLastError")
         import
      end function cudaGetLastError

      integer(c_int) function cudaGetMipmappedArrayLevel(levelArray, mipmappedArray, level) &
         bind(C, name="cudaGetMipmappedArrayLevel")
         import
         type(c_ptr), intent(out) :: levelArray
         type(c_ptr), value :: mipmappedArray
         integer(c_int), value :: level
      end function cudaGetMipmappedArrayLevel

      integer(c_int) function cudaGetSurfaceObjectResourceDesc(pResDesc, surfObject) &
         bind(C, name="cudaGetSurfaceObjectResourceDesc")
         import
         type(cudaResourceDesc), intent(inout) :: pResDesc
         integer(c_long_long), value :: surfObject
      end function cudaGetSurfaceObjectResourceDesc

      integer(c_int) function cudaGetSymbolAddress(devPtr, symbol) &
         bind(C, name="cudaGetSymbolAddress")
         import
         type(c_ptr), intent(out) :: devPtr
         type(c_ptr), value :: symbol
      end function cudaGetSymbolAddress

      integer(c_int) function cudaGetSymbolSize(size, symbol) &
         bind(C, name="cudaGetSymbolSize")
         import
         integer(c_size_t), intent(inout) :: size
         type(c_ptr), value :: symbol
      end function cudaGetSymbolSize

      integer(c_int) function cudaGetTextureObjectResourceDesc(pResDesc, texObject) &
         bind(C, name="cudaGetTextureObjectResourceDesc")
         import
         type(cudaResourceDesc), intent(inout) :: pResDesc
         integer(c_long_long), value :: texObject
      end function cudaGetTextureObjectResourceDesc

      integer(c_int) function cudaGetTextureObjectResourceViewDesc(pResViewDesc, texObject) &
         bind(C, name="cudaGetTextureObjectResourceViewDesc")
         import
         type(cudaResourceViewDesc), intent(inout) :: pResViewDesc
         integer(c_long_long), value :: texObject
      end function cudaGetTextureObjectResourceViewDesc

      integer(c_int) function cudaGetTextureObjectTextureDesc(pTexDesc, texObject) &
         bind(C, name="cudaGetTextureObjectTextureDesc")
         import
         type(cudaTextureDesc), intent(inout) :: pTexDesc
         integer(c_long_long), value :: texObject
      end function cudaGetTextureObjectTextureDesc

      integer(c_int) function cudaGraphAddChildGraphNode( &
         pGraphNode, graph, pDependencies, numDependencies, childGraph) &
         bind(C, name="cudaGraphAddChildGraphNode")
         import
         type(c_ptr), intent(out) :: pGraphNode
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: pDependencies
         integer(c_size_t), value :: numDependencies
         type(c_ptr), value :: childGraph
      end function cudaGraphAddChildGraphNode

      integer(c_int) function cudaGraphAddDependencies(graph, from, to, numDependencies) &
         bind(C, name="cudaGraphAddDependencies")
         import
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: from
         type(c_ptr), intent(out) :: to
         integer(c_size_t), value :: numDependencies
      end function cudaGraphAddDependencies

      integer(c_int) function cudaGraphAddDependencies_v2(graph, from, to, edgeData, numDependencies) &
         bind(C, name="cudaGraphAddDependencies_v2")
         import
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: from
         type(c_ptr), intent(out) :: to
         type(cudaGraphEdgeData_st), intent(in) :: edgeData
         integer(c_size_t), value :: numDependencies
      end function cudaGraphAddDependencies_v2

      integer(c_int) function cudaGraphAddEmptyNode(pGraphNode, graph, pDependencies, numDependencies) &
         bind(C, name="cudaGraphAddEmptyNode")
         import
         type(c_ptr), intent(out) :: pGraphNode
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: pDependencies
         integer(c_size_t), value :: numDependencies
      end function cudaGraphAddEmptyNode

      integer(c_int) function cudaGraphAddEventRecordNode(pGraphNode, graph, pDependencies, numDependencies, event) &
         bind(C, name="cudaGraphAddEventRecordNode")
         import
         type(c_ptr), intent(out) :: pGraphNode
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: pDependencies
         integer(c_size_t), value :: numDependencies
         type(c_ptr), value :: event
      end function cudaGraphAddEventRecordNode

      integer(c_int) function cudaGraphAddEventWaitNode(pGraphNode, graph, pDependencies, numDependencies, event) &
         bind(C, name="cudaGraphAddEventWaitNode")
         import
         type(c_ptr), intent(out) :: pGraphNode
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: pDependencies
         integer(c_size_t), value :: numDependencies
         type(c_ptr), value :: event
      end function cudaGraphAddEventWaitNode

      integer(c_int) function cudaGraphAddExternalSemaphoresSignalNode( &
         pGraphNode, graph, pDependencies, numDependencies, nodeParams) &
         bind(C, name="cudaGraphAddExternalSemaphoresSignalNode")
         import
         type(c_ptr), intent(out) :: pGraphNode
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: pDependencies
         integer(c_size_t), value :: numDependencies
         type(cudaExternalSemaphoreSignalNodeParams), intent(in) :: nodeParams
      end function cudaGraphAddExternalSemaphoresSignalNode

      integer(c_int) function cudaGraphAddExternalSemaphoresWaitNode( &
         pGraphNode, graph, pDependencies, numDependencies, nodeParams) &
         bind(C, name="cudaGraphAddExternalSemaphoresWaitNode")
         import
         type(c_ptr), intent(out) :: pGraphNode
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: pDependencies
         integer(c_size_t), value :: numDependencies
         type(cudaExternalSemaphoreWaitNodeParams), intent(in) :: nodeParams
      end function cudaGraphAddExternalSemaphoresWaitNode

      integer(c_int) function cudaGraphAddHostNode(pGraphNode, graph, pDependencies, numDependencies, pNodeParams) &
         bind(C, name="cudaGraphAddHostNode")
         import
         type(c_ptr), intent(out) :: pGraphNode
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: pDependencies
         integer(c_size_t), value :: numDependencies
         type(cudaHostNodeParams), intent(in) :: pNodeParams
      end function cudaGraphAddHostNode

      integer(c_int) function cudaGraphAddKernelNode(pGraphNode, graph, pDependencies, numDependencies, pNodeParams) &
         bind(C, name="cudaGraphAddKernelNode")
         import
         type(c_ptr), intent(out) :: pGraphNode
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: pDependencies
         integer(c_size_t), value :: numDependencies
         type(cudaKernelNodeParams), intent(in) :: pNodeParams
      end function cudaGraphAddKernelNode

      integer(c_int) function cudaGraphAddMemAllocNode( &
         pGraphNode, graph, pDependencies, numDependencies, nodeParams) &
         bind(C, name="cudaGraphAddMemAllocNode")
         import
         type(c_ptr), intent(out) :: pGraphNode
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: pDependencies
         integer(c_size_t), value :: numDependencies
         type(cudaMemAllocNodeParams), intent(inout) :: nodeParams
      end function cudaGraphAddMemAllocNode

      integer(c_int) function cudaGraphAddMemFreeNode(pGraphNode, graph, pDependencies, numDependencies, dptr) &
         bind(C, name="cudaGraphAddMemFreeNode")
         import
         type(c_ptr), intent(out) :: pGraphNode
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: pDependencies
         integer(c_size_t), value :: numDependencies
         type(c_ptr), value :: dptr
      end function cudaGraphAddMemFreeNode

      integer(c_int) function cudaGraphAddMemcpyNode(pGraphNode, graph, pDependencies, numDependencies, pCopyParams) &
         bind(C, name="cudaGraphAddMemcpyNode")
         import
         type(c_ptr), intent(out) :: pGraphNode
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: pDependencies
         integer(c_size_t), value :: numDependencies
         type(cudaMemcpy3DParms), intent(in) :: pCopyParams
      end function cudaGraphAddMemcpyNode

      integer(c_int) function cudaGraphAddMemcpyNode1D( &
         pGraphNode, graph, pDependencies, numDependencies, dst, src, count, kind) &
         bind(C, name="cudaGraphAddMemcpyNode1D")
         import
         type(c_ptr), intent(out) :: pGraphNode
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: pDependencies
         integer(c_size_t), value :: numDependencies
         type(c_ptr), value :: dst
         type(c_ptr), value :: src
         integer(c_size_t), value :: count
         integer(c_int), value :: kind
      end function cudaGraphAddMemcpyNode1D

      integer(c_int) function cudaGraphAddMemcpyNodeFromSymbol( &
         pGraphNode, graph, pDependencies, numDependencies, dst, symbol, count, offset, kind) &
         bind(C, name="cudaGraphAddMemcpyNodeFromSymbol")
         import
         type(c_ptr), intent(out) :: pGraphNode
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: pDependencies
         integer(c_size_t), value :: numDependencies
         type(c_ptr), value :: dst
         type(c_ptr), value :: symbol
         integer(c_size_t), value :: count
         integer(c_size_t), value :: offset
         integer(c_int), value :: kind
      end function cudaGraphAddMemcpyNodeFromSymbol

      integer(c_int) function cudaGraphAddMemcpyNodeToSymbol( &
         pGraphNode, graph, pDependencies, numDependencies, symbol, src, count, offset, kind) &
         bind(C, name="cudaGraphAddMemcpyNodeToSymbol")
         import
         type(c_ptr), intent(out) :: pGraphNode
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: pDependencies
         integer(c_size_t), value :: numDependencies
         type(c_ptr), value :: symbol
         type(c_ptr), value :: src
         integer(c_size_t), value :: count
         integer(c_size_t), value :: offset
         integer(c_int), value :: kind
      end function cudaGraphAddMemcpyNodeToSymbol

      integer(c_int) function cudaGraphAddMemsetNode( &
         pGraphNode, graph, pDependencies, numDependencies, pMemsetParams) &
         bind(C, name="cudaGraphAddMemsetNode")
         import
         type(c_ptr), intent(out) :: pGraphNode
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: pDependencies
         integer(c_size_t), value :: numDependencies
         type(cudaMemsetParams), intent(in) :: pMemsetParams
      end function cudaGraphAddMemsetNode

      integer(c_int) function cudaGraphAddNode(pGraphNode, graph, pDependencies, numDependencies, nodeParams) &
         bind(C, name="cudaGraphAddNode")
         import
         type(c_ptr), intent(out) :: pGraphNode
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: pDependencies
         integer(c_size_t), value :: numDependencies
         type(cudaGraphNodeParams), intent(inout) :: nodeParams
      end function cudaGraphAddNode

      integer(c_int) function cudaGraphAddNode_v2( &
         pGraphNode, graph, pDependencies, dependencyData, numDependencies, nodeParams) &
         bind(C, name="cudaGraphAddNode_v2")
         import
         type(c_ptr), intent(out) :: pGraphNode
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: pDependencies
         type(cudaGraphEdgeData_st), intent(in) :: dependencyData
         integer(c_size_t), value :: numDependencies
         type(cudaGraphNodeParams), intent(inout) :: nodeParams
      end function cudaGraphAddNode_v2

      integer(c_int) function cudaGraphChildGraphNodeGetGraph(node, pGraph) &
         bind(C, name="cudaGraphChildGraphNodeGetGraph")
         import
         type(c_ptr), value :: node
         type(c_ptr), intent(out) :: pGraph
      end function cudaGraphChildGraphNodeGetGraph

      integer(c_int) function cudaGraphClone(pGraphClone, originalGraph) &
         bind(C, name="cudaGraphClone")
         import
         type(c_ptr), intent(out) :: pGraphClone
         type(c_ptr), value :: originalGraph
      end function cudaGraphClone

      integer(c_int) function cudaGraphConditionalHandleCreate(pHandle_out, graph, defaultLaunchValue, flags) &
         bind(C, name="cudaGraphConditionalHandleCreate")
         import
         integer(c_long_long), intent(inout) :: pHandle_out
         type(c_ptr), value :: graph
         integer(c_int), value :: defaultLaunchValue
         integer(c_int), value :: flags
      end function cudaGraphConditionalHandleCreate

      integer(c_int) function cudaGraphCreate(pGraph, flags) &
         bind(C, name="cudaGraphCreate")
         import
         type(c_ptr), intent(out) :: pGraph
         integer(c_int), value :: flags
      end function cudaGraphCreate

      integer(c_int) function cudaGraphDebugDotPrint(graph, path, flags) &
         bind(C, name="cudaGraphDebugDotPrint")
         import
         type(c_ptr), value :: graph
         character(kind=c_char), dimension(*), intent(in) :: path
         integer(c_int), value :: flags
      end function cudaGraphDebugDotPrint

      integer(c_int) function cudaGraphDestroy(graph) &
         bind(C, name="cudaGraphDestroy")
         import
         type(c_ptr), value :: graph
      end function cudaGraphDestroy

      integer(c_int) function cudaGraphDestroyNode(node) &
         bind(C, name="cudaGraphDestroyNode")
         import
         type(c_ptr), value :: node
      end function cudaGraphDestroyNode

      integer(c_int) function cudaGraphEventRecordNodeGetEvent(node, event_out) &
         bind(C, name="cudaGraphEventRecordNodeGetEvent")
         import
         type(c_ptr), value :: node
         type(c_ptr), intent(out) :: event_out
      end function cudaGraphEventRecordNodeGetEvent

      integer(c_int) function cudaGraphEventRecordNodeSetEvent(node, event) &
         bind(C, name="cudaGraphEventRecordNodeSetEvent")
         import
         type(c_ptr), value :: node
         type(c_ptr), value :: event
      end function cudaGraphEventRecordNodeSetEvent

      integer(c_int) function cudaGraphEventWaitNodeGetEvent(node, event_out) &
         bind(C, name="cudaGraphEventWaitNodeGetEvent")
         import
         type(c_ptr), value :: node
         type(c_ptr), intent(out) :: event_out
      end function cudaGraphEventWaitNodeGetEvent

      integer(c_int) function cudaGraphEventWaitNodeSetEvent(node, event) &
         bind(C, name="cudaGraphEventWaitNodeSetEvent")
         import
         type(c_ptr), value :: node
         type(c_ptr), value :: event
      end function cudaGraphEventWaitNodeSetEvent

      integer(c_int) function cudaGraphExecChildGraphNodeSetParams(hGraphExec, node, childGraph) &
         bind(C, name="cudaGraphExecChildGraphNodeSetParams")
         import
         type(c_ptr), value :: hGraphExec
         type(c_ptr), value :: node
         type(c_ptr), value :: childGraph
      end function cudaGraphExecChildGraphNodeSetParams

      integer(c_int) function cudaGraphExecDestroy(graphExec) &
         bind(C, name="cudaGraphExecDestroy")
         import
         type(c_ptr), value :: graphExec
      end function cudaGraphExecDestroy

      integer(c_int) function cudaGraphExecEventRecordNodeSetEvent(hGraphExec, hNode, event) &
         bind(C, name="cudaGraphExecEventRecordNodeSetEvent")
         import
         type(c_ptr), value :: hGraphExec
         type(c_ptr), value :: hNode
         type(c_ptr), value :: event
      end function cudaGraphExecEventRecordNodeSetEvent

      integer(c_int) function cudaGraphExecEventWaitNodeSetEvent(hGraphExec, hNode, event) &
         bind(C, name="cudaGraphExecEventWaitNodeSetEvent")
         import
         type(c_ptr), value :: hGraphExec
         type(c_ptr), value :: hNode
         type(c_ptr), value :: event
      end function cudaGraphExecEventWaitNodeSetEvent

      integer(c_int) function cudaGraphExecExternalSemaphoresSignalNodeSetParams(hGraphExec, hNode, nodeParams) &
         bind(C, name="cudaGraphExecExternalSemaphoresSignalNodeSetParams")
         import
         type(c_ptr), value :: hGraphExec
         type(c_ptr), value :: hNode
         type(cudaExternalSemaphoreSignalNodeParams), intent(in) :: nodeParams
      end function cudaGraphExecExternalSemaphoresSignalNodeSetParams

      integer(c_int) function cudaGraphExecExternalSemaphoresWaitNodeSetParams(hGraphExec, hNode, nodeParams) &
         bind(C, name="cudaGraphExecExternalSemaphoresWaitNodeSetParams")
         import
         type(c_ptr), value :: hGraphExec
         type(c_ptr), value :: hNode
         type(cudaExternalSemaphoreWaitNodeParams), intent(in) :: nodeParams
      end function cudaGraphExecExternalSemaphoresWaitNodeSetParams

      integer(c_int) function cudaGraphExecGetFlags(graphExec, flags) &
         bind(C, name="cudaGraphExecGetFlags")
         import
         type(c_ptr), value :: graphExec
         integer(c_long_long), intent(inout) :: flags
      end function cudaGraphExecGetFlags

      integer(c_int) function cudaGraphExecHostNodeSetParams(hGraphExec, node, pNodeParams) &
         bind(C, name="cudaGraphExecHostNodeSetParams")
         import
         type(c_ptr), value :: hGraphExec
         type(c_ptr), value :: node
         type(cudaHostNodeParams), intent(in) :: pNodeParams
      end function cudaGraphExecHostNodeSetParams

      integer(c_int) function cudaGraphExecKernelNodeSetParams(hGraphExec, node, pNodeParams) &
         bind(C, name="cudaGraphExecKernelNodeSetParams")
         import
         type(c_ptr), value :: hGraphExec
         type(c_ptr), value :: node
         type(cudaKernelNodeParams), intent(in) :: pNodeParams
      end function cudaGraphExecKernelNodeSetParams

      integer(c_int) function cudaGraphExecMemcpyNodeSetParams(hGraphExec, node, pNodeParams) &
         bind(C, name="cudaGraphExecMemcpyNodeSetParams")
         import
         type(c_ptr), value :: hGraphExec
         type(c_ptr), value :: node
         type(cudaMemcpy3DParms), intent(in) :: pNodeParams
      end function cudaGraphExecMemcpyNodeSetParams

      integer(c_int) function cudaGraphExecMemcpyNodeSetParams1D(hGraphExec, node, dst, src, count, kind) &
         bind(C, name="cudaGraphExecMemcpyNodeSetParams1D")
         import
         type(c_ptr), value :: hGraphExec
         type(c_ptr), value :: node
         type(c_ptr), value :: dst
         type(c_ptr), value :: src
         integer(c_size_t), value :: count
         integer(c_int), value :: kind
      end function cudaGraphExecMemcpyNodeSetParams1D

      integer(c_int) function cudaGraphExecMemcpyNodeSetParamsFromSymbol( &
         hGraphExec, node, dst, symbol, count, offset, kind) &
         bind(C, name="cudaGraphExecMemcpyNodeSetParamsFromSymbol")
         import
         type(c_ptr), value :: hGraphExec
         type(c_ptr), value :: node
         type(c_ptr), value :: dst
         type(c_ptr), value :: symbol
         integer(c_size_t), value :: count
         integer(c_size_t), value :: offset
         integer(c_int), value :: kind
      end function cudaGraphExecMemcpyNodeSetParamsFromSymbol

      integer(c_int) function cudaGraphExecMemcpyNodeSetParamsToSymbol( &
         hGraphExec, node, symbol, src, count, offset, kind) &
         bind(C, name="cudaGraphExecMemcpyNodeSetParamsToSymbol")
         import
         type(c_ptr), value :: hGraphExec
         type(c_ptr), value :: node
         type(c_ptr), value :: symbol
         type(c_ptr), value :: src
         integer(c_size_t), value :: count
         integer(c_size_t), value :: offset
         integer(c_int), value :: kind
      end function cudaGraphExecMemcpyNodeSetParamsToSymbol

      integer(c_int) function cudaGraphExecMemsetNodeSetParams(hGraphExec, node, pNodeParams) &
         bind(C, name="cudaGraphExecMemsetNodeSetParams")
         import
         type(c_ptr), value :: hGraphExec
         type(c_ptr), value :: node
         type(cudaMemsetParams), intent(in) :: pNodeParams
      end function cudaGraphExecMemsetNodeSetParams

      integer(c_int) function cudaGraphExecNodeSetParams(graphExec, node, nodeParams) &
         bind(C, name="cudaGraphExecNodeSetParams")
         import
         type(c_ptr), value :: graphExec
         type(c_ptr), value :: node
         type(cudaGraphNodeParams), intent(inout) :: nodeParams
      end function cudaGraphExecNodeSetParams

      integer(c_int) function cudaGraphExecUpdate(hGraphExec, hGraph, resultInfo) &
         bind(C, name="cudaGraphExecUpdate")
         import
         type(c_ptr), value :: hGraphExec
         type(c_ptr), value :: hGraph
         type(cudaGraphExecUpdateResultInfo_st), intent(inout) :: resultInfo
      end function cudaGraphExecUpdate

      integer(c_int) function cudaGraphExternalSemaphoresSignalNodeGetParams(hNode, params_out) &
         bind(C, name="cudaGraphExternalSemaphoresSignalNodeGetParams")
         import
         type(c_ptr), value :: hNode
         type(cudaExternalSemaphoreSignalNodeParams), intent(inout) :: params_out
      end function cudaGraphExternalSemaphoresSignalNodeGetParams

      integer(c_int) function cudaGraphExternalSemaphoresSignalNodeSetParams(hNode, nodeParams) &
         bind(C, name="cudaGraphExternalSemaphoresSignalNodeSetParams")
         import
         type(c_ptr), value :: hNode
         type(cudaExternalSemaphoreSignalNodeParams), intent(in) :: nodeParams
      end function cudaGraphExternalSemaphoresSignalNodeSetParams

      integer(c_int) function cudaGraphExternalSemaphoresWaitNodeGetParams(hNode, params_out) &
         bind(C, name="cudaGraphExternalSemaphoresWaitNodeGetParams")
         import
         type(c_ptr), value :: hNode
         type(cudaExternalSemaphoreWaitNodeParams), intent(inout) :: params_out
      end function cudaGraphExternalSemaphoresWaitNodeGetParams

      integer(c_int) function cudaGraphExternalSemaphoresWaitNodeSetParams(hNode, nodeParams) &
         bind(C, name="cudaGraphExternalSemaphoresWaitNodeSetParams")
         import
         type(c_ptr), value :: hNode
         type(cudaExternalSemaphoreWaitNodeParams), intent(in) :: nodeParams
      end function cudaGraphExternalSemaphoresWaitNodeSetParams

      integer(c_int) function cudaGraphGetEdges(graph, from, to, numEdges) &
         bind(C, name="cudaGraphGetEdges")
         import
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: from
         type(c_ptr), intent(out) :: to
         integer(c_size_t), intent(inout) :: numEdges
      end function cudaGraphGetEdges

      integer(c_int) function cudaGraphGetEdges_v2(graph, from, to, edgeData, numEdges) &
         bind(C, name="cudaGraphGetEdges_v2")
         import
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: from
         type(c_ptr), intent(out) :: to
         type(cudaGraphEdgeData_st), intent(inout) :: edgeData
         integer(c_size_t), intent(inout) :: numEdges
      end function cudaGraphGetEdges_v2

      integer(c_int) function cudaGraphGetNodes(graph, nodes, numNodes) &
         bind(C, name="cudaGraphGetNodes")
         import
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: nodes
         integer(c_size_t), intent(inout) :: numNodes
      end function cudaGraphGetNodes

      integer(c_int) function cudaGraphGetRootNodes(graph, pRootNodes, pNumRootNodes) &
         bind(C, name="cudaGraphGetRootNodes")
         import
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: pRootNodes
         integer(c_size_t), intent(inout) :: pNumRootNodes
      end function cudaGraphGetRootNodes

      integer(c_int) function cudaGraphHostNodeGetParams(node, pNodeParams) &
         bind(C, name="cudaGraphHostNodeGetParams")
         import
         type(c_ptr), value :: node
         type(cudaHostNodeParams), intent(inout) :: pNodeParams
      end function cudaGraphHostNodeGetParams

      integer(c_int) function cudaGraphHostNodeSetParams(node, pNodeParams) &
         bind(C, name="cudaGraphHostNodeSetParams")
         import
         type(c_ptr), value :: node
         type(cudaHostNodeParams), intent(in) :: pNodeParams
      end function cudaGraphHostNodeSetParams

      integer(c_int) function cudaGraphInstantiate(pGraphExec, graph, flags) &
         bind(C, name="cudaGraphInstantiate")
         import
         type(c_ptr), intent(out) :: pGraphExec
         type(c_ptr), value :: graph
         integer(c_long_long), value :: flags
      end function cudaGraphInstantiate

      integer(c_int) function cudaGraphInstantiateWithFlags(pGraphExec, graph, flags) &
         bind(C, name="cudaGraphInstantiateWithFlags")
         import
         type(c_ptr), intent(out) :: pGraphExec
         type(c_ptr), value :: graph
         integer(c_long_long), value :: flags
      end function cudaGraphInstantiateWithFlags

      integer(c_int) function cudaGraphInstantiateWithParams(pGraphExec, graph, instantiateParams) &
         bind(C, name="cudaGraphInstantiateWithParams")
         import
         type(c_ptr), intent(out) :: pGraphExec
         type(c_ptr), value :: graph
         type(cudaGraphInstantiateParams_st), intent(inout) :: instantiateParams
      end function cudaGraphInstantiateWithParams

      integer(c_int) function cudaGraphKernelNodeCopyAttributes(hSrc, hDst) &
         bind(C, name="cudaGraphKernelNodeCopyAttributes")
         import
         type(c_ptr), value :: hSrc
         type(c_ptr), value :: hDst
      end function cudaGraphKernelNodeCopyAttributes

      integer(c_int) function cudaGraphKernelNodeGetAttribute(hNode, attr, value_out) &
         bind(C, name="cudaGraphKernelNodeGetAttribute")
         import
         type(c_ptr), value :: hNode
         integer(c_int), value :: attr
         type(cudaLaunchAttributeValue), intent(inout) :: value_out
      end function cudaGraphKernelNodeGetAttribute

      integer(c_int) function cudaGraphKernelNodeGetParams(node, pNodeParams) &
         bind(C, name="cudaGraphKernelNodeGetParams")
         import
         type(c_ptr), value :: node
         type(cudaKernelNodeParams), intent(inout) :: pNodeParams
      end function cudaGraphKernelNodeGetParams

      integer(c_int) function cudaGraphKernelNodeSetAttribute(hNode, attr, value) &
         bind(C, name="cudaGraphKernelNodeSetAttribute")
         import
         type(c_ptr), value :: hNode
         integer(c_int), value :: attr
         type(cudaLaunchAttributeValue), intent(in) :: value
      end function cudaGraphKernelNodeSetAttribute

      integer(c_int) function cudaGraphKernelNodeSetParams(node, pNodeParams) &
         bind(C, name="cudaGraphKernelNodeSetParams")
         import
         type(c_ptr), value :: node
         type(cudaKernelNodeParams), intent(in) :: pNodeParams
      end function cudaGraphKernelNodeSetParams

      integer(c_int) function cudaGraphLaunch(graphExec, stream) &
         bind(C, name="cudaGraphLaunch")
         import
         type(c_ptr), value :: graphExec
         type(c_ptr), value :: stream
      end function cudaGraphLaunch

      integer(c_int) function cudaGraphMemAllocNodeGetParams(node, params_out) &
         bind(C, name="cudaGraphMemAllocNodeGetParams")
         import
         type(c_ptr), value :: node
         type(cudaMemAllocNodeParams), intent(inout) :: params_out
      end function cudaGraphMemAllocNodeGetParams

      integer(c_int) function cudaGraphMemFreeNodeGetParams(node, dptr_out) &
         bind(C, name="cudaGraphMemFreeNodeGetParams")
         import
         type(c_ptr), value :: node
         type(c_ptr), value :: dptr_out
      end function cudaGraphMemFreeNodeGetParams

      integer(c_int) function cudaGraphMemcpyNodeGetParams(node, pNodeParams) &
         bind(C, name="cudaGraphMemcpyNodeGetParams")
         import
         type(c_ptr), value :: node
         type(cudaMemcpy3DParms), intent(inout) :: pNodeParams
      end function cudaGraphMemcpyNodeGetParams

      integer(c_int) function cudaGraphMemcpyNodeSetParams(node, pNodeParams) &
         bind(C, name="cudaGraphMemcpyNodeSetParams")
         import
         type(c_ptr), value :: node
         type(cudaMemcpy3DParms), intent(in) :: pNodeParams
      end function cudaGraphMemcpyNodeSetParams

      integer(c_int) function cudaGraphMemcpyNodeSetParams1D(node, dst, src, count, kind) &
         bind(C, name="cudaGraphMemcpyNodeSetParams1D")
         import
         type(c_ptr), value :: node
         type(c_ptr), value :: dst
         type(c_ptr), value :: src
         integer(c_size_t), value :: count
         integer(c_int), value :: kind
      end function cudaGraphMemcpyNodeSetParams1D

      integer(c_int) function cudaGraphMemcpyNodeSetParamsFromSymbol(node, dst, symbol, count, offset, kind) &
         bind(C, name="cudaGraphMemcpyNodeSetParamsFromSymbol")
         import
         type(c_ptr), value :: node
         type(c_ptr), value :: dst
         type(c_ptr), value :: symbol
         integer(c_size_t), value :: count
         integer(c_size_t), value :: offset
         integer(c_int), value :: kind
      end function cudaGraphMemcpyNodeSetParamsFromSymbol

      integer(c_int) function cudaGraphMemcpyNodeSetParamsToSymbol(node, symbol, src, count, offset, kind) &
         bind(C, name="cudaGraphMemcpyNodeSetParamsToSymbol")
         import
         type(c_ptr), value :: node
         type(c_ptr), value :: symbol
         type(c_ptr), value :: src
         integer(c_size_t), value :: count
         integer(c_size_t), value :: offset
         integer(c_int), value :: kind
      end function cudaGraphMemcpyNodeSetParamsToSymbol

      integer(c_int) function cudaGraphMemsetNodeGetParams(node, pNodeParams) &
         bind(C, name="cudaGraphMemsetNodeGetParams")
         import
         type(c_ptr), value :: node
         type(cudaMemsetParams), intent(inout) :: pNodeParams
      end function cudaGraphMemsetNodeGetParams

      integer(c_int) function cudaGraphMemsetNodeSetParams(node, pNodeParams) &
         bind(C, name="cudaGraphMemsetNodeSetParams")
         import
         type(c_ptr), value :: node
         type(cudaMemsetParams), intent(in) :: pNodeParams
      end function cudaGraphMemsetNodeSetParams

      integer(c_int) function cudaGraphNodeFindInClone(pNode, originalNode, clonedGraph) &
         bind(C, name="cudaGraphNodeFindInClone")
         import
         type(c_ptr), intent(out) :: pNode
         type(c_ptr), value :: originalNode
         type(c_ptr), value :: clonedGraph
      end function cudaGraphNodeFindInClone

      integer(c_int) function cudaGraphNodeGetDependencies(node, pDependencies, pNumDependencies) &
         bind(C, name="cudaGraphNodeGetDependencies")
         import
         type(c_ptr), value :: node
         type(c_ptr), intent(out) :: pDependencies
         integer(c_size_t), intent(inout) :: pNumDependencies
      end function cudaGraphNodeGetDependencies

      integer(c_int) function cudaGraphNodeGetDependencies_v2(node, pDependencies, edgeData, pNumDependencies) &
         bind(C, name="cudaGraphNodeGetDependencies_v2")
         import
         type(c_ptr), value :: node
         type(c_ptr), intent(out) :: pDependencies
         type(cudaGraphEdgeData_st), intent(inout) :: edgeData
         integer(c_size_t), intent(inout) :: pNumDependencies
      end function cudaGraphNodeGetDependencies_v2

      integer(c_int) function cudaGraphNodeGetDependentNodes(node, pDependentNodes, pNumDependentNodes) &
         bind(C, name="cudaGraphNodeGetDependentNodes")
         import
         type(c_ptr), value :: node
         type(c_ptr), intent(out) :: pDependentNodes
         integer(c_size_t), intent(inout) :: pNumDependentNodes
      end function cudaGraphNodeGetDependentNodes

      integer(c_int) function cudaGraphNodeGetDependentNodes_v2(node, pDependentNodes, edgeData, pNumDependentNodes) &
         bind(C, name="cudaGraphNodeGetDependentNodes_v2")
         import
         type(c_ptr), value :: node
         type(c_ptr), intent(out) :: pDependentNodes
         type(cudaGraphEdgeData_st), intent(inout) :: edgeData
         integer(c_size_t), intent(inout) :: pNumDependentNodes
      end function cudaGraphNodeGetDependentNodes_v2

      integer(c_int) function cudaGraphNodeGetEnabled(hGraphExec, hNode, isEnabled) &
         bind(C, name="cudaGraphNodeGetEnabled")
         import
         type(c_ptr), value :: hGraphExec
         type(c_ptr), value :: hNode
         integer(c_int), intent(inout) :: isEnabled
      end function cudaGraphNodeGetEnabled

      integer(c_int) function cudaGraphNodeGetType(node, pType) &
         bind(C, name="cudaGraphNodeGetType")
         import
         type(c_ptr), value :: node
         integer(c_int), intent(out) :: pType
      end function cudaGraphNodeGetType

      integer(c_int) function cudaGraphNodeSetEnabled(hGraphExec, hNode, isEnabled) &
         bind(C, name="cudaGraphNodeSetEnabled")
         import
         type(c_ptr), value :: hGraphExec
         type(c_ptr), value :: hNode
         integer(c_int), value :: isEnabled
      end function cudaGraphNodeSetEnabled

      integer(c_int) function cudaGraphNodeSetParams(node, nodeParams) &
         bind(C, name="cudaGraphNodeSetParams")
         import
         type(c_ptr), value :: node
         type(cudaGraphNodeParams), intent(inout) :: nodeParams
      end function cudaGraphNodeSetParams

      integer(c_int) function cudaGraphReleaseUserObject(graph, object, count) &
         bind(C, name="cudaGraphReleaseUserObject")
         import
         type(c_ptr), value :: graph
         type(c_ptr), value :: object
         integer(c_int), value :: count
      end function cudaGraphReleaseUserObject

      integer(c_int) function cudaGraphRemoveDependencies(graph, from, to, numDependencies) &
         bind(C, name="cudaGraphRemoveDependencies")
         import
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: from
         type(c_ptr), intent(out) :: to
         integer(c_size_t), value :: numDependencies
      end function cudaGraphRemoveDependencies

      integer(c_int) function cudaGraphRemoveDependencies_v2(graph, from, to, edgeData, numDependencies) &
         bind(C, name="cudaGraphRemoveDependencies_v2")
         import
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: from
         type(c_ptr), intent(out) :: to
         type(cudaGraphEdgeData_st), intent(in) :: edgeData
         integer(c_size_t), value :: numDependencies
      end function cudaGraphRemoveDependencies_v2

      integer(c_int) function cudaGraphRetainUserObject(graph, object, count, flags) &
         bind(C, name="cudaGraphRetainUserObject")
         import
         type(c_ptr), value :: graph
         type(c_ptr), value :: object
         integer(c_int), value :: count
         integer(c_int), value :: flags
      end function cudaGraphRetainUserObject

      integer(c_int) function cudaGraphUpload(graphExec, stream) &
         bind(C, name="cudaGraphUpload")
         import
         type(c_ptr), value :: graphExec
         type(c_ptr), value :: stream
      end function cudaGraphUpload

      integer(c_int) function cudaGraphicsMapResources(count, resources, stream) &
         bind(C, name="cudaGraphicsMapResources")
         import
         integer(c_int), value :: count
         type(c_ptr), intent(out) :: resources
         type(c_ptr), value :: stream
      end function cudaGraphicsMapResources

      integer(c_int) function cudaGraphicsResourceGetMappedMipmappedArray(mipmappedArray, resource) &
         bind(C, name="cudaGraphicsResourceGetMappedMipmappedArray")
         import
         type(c_ptr), intent(out) :: mipmappedArray
         type(c_ptr), value :: resource
      end function cudaGraphicsResourceGetMappedMipmappedArray

      integer(c_int) function cudaGraphicsResourceGetMappedPointer(devPtr, size, resource) &
         bind(C, name="cudaGraphicsResourceGetMappedPointer")
         import
         type(c_ptr), intent(out) :: devPtr
         integer(c_size_t), intent(inout) :: size
         type(c_ptr), value :: resource
      end function cudaGraphicsResourceGetMappedPointer

      integer(c_int) function cudaGraphicsResourceSetMapFlags(resource, flags) &
         bind(C, name="cudaGraphicsResourceSetMapFlags")
         import
         type(c_ptr), value :: resource
         integer(c_int), value :: flags
      end function cudaGraphicsResourceSetMapFlags

      integer(c_int) function cudaGraphicsSubResourceGetMappedArray(array, resource, arrayIndex, mipLevel) &
         bind(C, name="cudaGraphicsSubResourceGetMappedArray")
         import
         type(c_ptr), intent(out) :: array
         type(c_ptr), value :: resource
         integer(c_int), value :: arrayIndex
         integer(c_int), value :: mipLevel
      end function cudaGraphicsSubResourceGetMappedArray

      integer(c_int) function cudaGraphicsUnmapResources(count, resources, stream) &
         bind(C, name="cudaGraphicsUnmapResources")
         import
         integer(c_int), value :: count
         type(c_ptr), intent(out) :: resources
         type(c_ptr), value :: stream
      end function cudaGraphicsUnmapResources

      integer(c_int) function cudaGraphicsUnregisterResource(resource) &
         bind(C, name="cudaGraphicsUnregisterResource")
         import
         type(c_ptr), value :: resource
      end function cudaGraphicsUnregisterResource

      integer(c_int) function cudaHostAlloc(pHost, size, flags) &
         bind(C, name="cudaHostAlloc")
         import
         type(c_ptr), intent(out) :: pHost
         integer(c_size_t), value :: size
         integer(c_int), value :: flags
      end function cudaHostAlloc

      integer(c_int) function cudaHostGetDevicePointer(pDevice, pHost, flags) &
         bind(C, name="cudaHostGetDevicePointer")
         import
         type(c_ptr), intent(out) :: pDevice
         type(c_ptr), value :: pHost
         integer(c_int), value :: flags
      end function cudaHostGetDevicePointer

      integer(c_int) function cudaHostGetFlags(pFlags, pHost) &
         bind(C, name="cudaHostGetFlags")
         import
         integer(c_int), intent(inout) :: pFlags
         type(c_ptr), value :: pHost
      end function cudaHostGetFlags

      integer(c_int) function cudaHostRegister(ptr, size, flags) &
         bind(C, name="cudaHostRegister")
         import
         type(c_ptr), value :: ptr
         integer(c_size_t), value :: size
         integer(c_int), value :: flags
      end function cudaHostRegister

      integer(c_int) function cudaHostUnregister(ptr) &
         bind(C, name="cudaHostUnregister")
         import
         type(c_ptr), value :: ptr
      end function cudaHostUnregister

      integer(c_int) function cudaImportExternalMemory(extMem_out, memHandleDesc) &
         bind(C, name="cudaImportExternalMemory")
         import
         type(c_ptr), intent(out) :: extMem_out
         type(cudaExternalMemoryHandleDesc), intent(in) :: memHandleDesc
      end function cudaImportExternalMemory

      integer(c_int) function cudaImportExternalSemaphore(extSem_out, semHandleDesc) &
         bind(C, name="cudaImportExternalSemaphore")
         import
         type(c_ptr), intent(out) :: extSem_out
         type(cudaExternalSemaphoreHandleDesc), intent(in) :: semHandleDesc
      end function cudaImportExternalSemaphore

      integer(c_int) function cudaInitDevice(device, deviceFlags, flags) &
         bind(C, name="cudaInitDevice")
         import
         integer(c_int), value :: device
         integer(c_int), value :: deviceFlags
         integer(c_int), value :: flags
      end function cudaInitDevice

      integer(c_int) function cudaIpcCloseMemHandle(devPtr) &
         bind(C, name="cudaIpcCloseMemHandle")
         import
         type(c_ptr), value :: devPtr
      end function cudaIpcCloseMemHandle

      integer(c_int) function cudaIpcGetEventHandle(handle, event) &
         bind(C, name="cudaIpcGetEventHandle")
         import
         type(cudaIpcEventHandle_st), intent(inout) :: handle
         type(c_ptr), value :: event
      end function cudaIpcGetEventHandle

      integer(c_int) function cudaIpcGetMemHandle(handle, devPtr) &
         bind(C, name="cudaIpcGetMemHandle")
         import
         type(cudaIpcMemHandle_st), intent(inout) :: handle
         type(c_ptr), value :: devPtr
      end function cudaIpcGetMemHandle

      integer(c_int) function cudaIpcOpenEventHandle(event, handle) &
         bind(C, name="cudaIpcOpenEventHandle")
         import
         type(c_ptr), intent(out) :: event
         type(cudaIpcEventHandle_st), value :: handle
      end function cudaIpcOpenEventHandle

      integer(c_int) function cudaIpcOpenMemHandle(devPtr, handle, flags) &
         bind(C, name="cudaIpcOpenMemHandle")
         import
         type(c_ptr), intent(out) :: devPtr
         type(cudaIpcMemHandle_st), value :: handle
         integer(c_int), value :: flags
      end function cudaIpcOpenMemHandle

      integer(c_int) function cudaKernelSetAttributeForDevice(kernel, attr, value, device) &
         bind(C, name="cudaKernelSetAttributeForDevice")
         import
         type(c_ptr), value :: kernel
         integer(c_int), value :: attr
         integer(c_int), value :: value
         integer(c_int), value :: device
      end function cudaKernelSetAttributeForDevice

      integer(c_int) function cudaLaunchCooperativeKernel(func, gridDim, blockDim, args, sharedMem, stream) &
         bind(C, name="cudaLaunchCooperativeKernel")
         import
         type(c_ptr), value :: func
         type(dim3), value :: gridDim
         type(dim3), value :: blockDim
         type(c_ptr), dimension(*), intent(in) :: args
         integer(c_size_t), value :: sharedMem
         type(c_ptr), value :: stream
      end function cudaLaunchCooperativeKernel

      integer(c_int) function cudaLaunchCooperativeKernelMultiDevice(launchParamsList, numDevices, flags) &
         bind(C, name="cudaLaunchCooperativeKernelMultiDevice")
         import
         type(cudaLaunchParams), intent(inout) :: launchParamsList
         integer(c_int), value :: numDevices
         integer(c_int), value :: flags
      end function cudaLaunchCooperativeKernelMultiDevice

      integer(c_int) function cudaLaunchHostFunc(stream, fn, userData) &
         bind(C, name="cudaLaunchHostFunc")
         import
         type(c_ptr), value :: stream
         type(c_funptr), value :: fn
         type(c_ptr), value :: userData
      end function cudaLaunchHostFunc

      integer(c_int) function cudaLaunchKernel(func, gridDim, blockDim, args, sharedMem, stream) &
         bind(C, name="cudaLaunchKernel")
         import
         type(c_ptr), value :: func
         type(dim3), value :: gridDim
         type(dim3), value :: blockDim
         type(c_ptr), dimension(*), intent(in) :: args
         integer(c_size_t), value :: sharedMem
         type(c_ptr), value :: stream
      end function cudaLaunchKernel

      integer(c_int) function cudaLaunchKernelExC(config, func, args) &
         bind(C, name="cudaLaunchKernelExC")
         import
         type(cudaLaunchConfig_st), intent(in) :: config
         type(c_ptr), value :: func
         type(c_ptr), dimension(*), intent(in) :: args
      end function cudaLaunchKernelExC

      integer(c_int) function cudaLibraryEnumerateKernels(kernels, numKernels, lib) &
         bind(C, name="cudaLibraryEnumerateKernels")
         import
         type(c_ptr), intent(out) :: kernels
         integer(c_int), value :: numKernels
         type(c_ptr), value :: lib
      end function cudaLibraryEnumerateKernels

      integer(c_int) function cudaLibraryGetGlobal(dptr, bytes, library, name) &
         bind(C, name="cudaLibraryGetGlobal")
         import
         type(c_ptr), intent(out) :: dptr
         integer(c_size_t), intent(inout) :: bytes
         type(c_ptr), value :: library
         character(kind=c_char), dimension(*), intent(in) :: name
      end function cudaLibraryGetGlobal

      integer(c_int) function cudaLibraryGetKernel(pKernel, library, name) &
         bind(C, name="cudaLibraryGetKernel")
         import
         type(c_ptr), intent(out) :: pKernel
         type(c_ptr), value :: library
         character(kind=c_char), dimension(*), intent(in) :: name
      end function cudaLibraryGetKernel

      integer(c_int) function cudaLibraryGetKernelCount(count, lib) &
         bind(C, name="cudaLibraryGetKernelCount")
         import
         integer(c_int), intent(inout) :: count
         type(c_ptr), value :: lib
      end function cudaLibraryGetKernelCount

      integer(c_int) function cudaLibraryGetManaged(dptr, bytes, library, name) &
         bind(C, name="cudaLibraryGetManaged")
         import
         type(c_ptr), intent(out) :: dptr
         integer(c_size_t), intent(inout) :: bytes
         type(c_ptr), value :: library
         character(kind=c_char), dimension(*), intent(in) :: name
      end function cudaLibraryGetManaged

      integer(c_int) function cudaLibraryGetUnifiedFunction(fptr, library, symbol) &
         bind(C, name="cudaLibraryGetUnifiedFunction")
         import
         type(c_ptr), intent(out) :: fptr
         type(c_ptr), value :: library
         character(kind=c_char), dimension(*), intent(in) :: symbol
      end function cudaLibraryGetUnifiedFunction

      integer(c_int) function cudaLibraryLoadData( &
         library, code, jitOptions, jitOptionsValues, numJitOptions, libraryOptions, libraryOptionValues, &
         numLibraryOptions) &
         bind(C, name="cudaLibraryLoadData")
         import
         type(c_ptr), intent(out) :: library
         type(c_ptr), value :: code
         integer(c_int), intent(out) :: jitOptions
         type(c_ptr), intent(out) :: jitOptionsValues
         integer(c_int), value :: numJitOptions
         integer(c_int), intent(out) :: libraryOptions
         type(c_ptr), intent(out) :: libraryOptionValues
         integer(c_int), value :: numLibraryOptions
      end function cudaLibraryLoadData

      integer(c_int) function cudaLibraryLoadFromFile( &
         library, fileName, jitOptions, jitOptionsValues, numJitOptions, libraryOptions, libraryOptionValues, &
         numLibraryOptions) &
         bind(C, name="cudaLibraryLoadFromFile")
         import
         type(c_ptr), intent(out) :: library
         character(kind=c_char), dimension(*), intent(in) :: fileName
         integer(c_int), intent(out) :: jitOptions
         type(c_ptr), intent(out) :: jitOptionsValues
         integer(c_int), value :: numJitOptions
         integer(c_int), intent(out) :: libraryOptions
         type(c_ptr), intent(out) :: libraryOptionValues
         integer(c_int), value :: numLibraryOptions
      end function cudaLibraryLoadFromFile

      integer(c_int) function cudaLibraryUnload(library) &
         bind(C, name="cudaLibraryUnload")
         import
         type(c_ptr), value :: library
      end function cudaLibraryUnload

      integer(c_int) function cudaMalloc(devPtr, size) &
         bind(C, name="cudaMalloc")
         import
         type(c_ptr), intent(out) :: devPtr
         integer(c_size_t), value :: size
      end function cudaMalloc

      integer(c_int) function cudaMalloc3D(pitchedDevPtr, extent) &
         bind(C, name="cudaMalloc3D")
         import
         type(cudaPitchedPtr), intent(inout) :: pitchedDevPtr
         type(cudaExtent), value :: extent
      end function cudaMalloc3D

      integer(c_int) function cudaMalloc3DArray(array, desc, extent, flags) &
         bind(C, name="cudaMalloc3DArray")
         import
         type(c_ptr), intent(out) :: array
         type(cudaChannelFormatDesc), intent(in) :: desc
         type(cudaExtent), value :: extent
         integer(c_int), value :: flags
      end function cudaMalloc3DArray

      integer(c_int) function cudaMallocArray(array, desc, width, height, flags) &
         bind(C, name="cudaMallocArray")
         import
         type(c_ptr), intent(out) :: array
         type(cudaChannelFormatDesc), intent(in) :: desc
         integer(c_size_t), value :: width
         integer(c_size_t), value :: height
         integer(c_int), value :: flags
      end function cudaMallocArray

      integer(c_int) function cudaMallocAsync(devPtr, size, hStream) &
         bind(C, name="cudaMallocAsync")
         import
         type(c_ptr), intent(out) :: devPtr
         integer(c_size_t), value :: size
         type(c_ptr), value :: hStream
      end function cudaMallocAsync

      integer(c_int) function cudaMallocFromPoolAsync(ptr, size, memPool, stream) &
         bind(C, name="cudaMallocFromPoolAsync")
         import
         type(c_ptr), intent(out) :: ptr
         integer(c_size_t), value :: size
         type(c_ptr), value :: memPool
         type(c_ptr), value :: stream
      end function cudaMallocFromPoolAsync

      integer(c_int) function cudaMallocHost(ptr, size) &
         bind(C, name="cudaMallocHost")
         import
         type(c_ptr), intent(out) :: ptr
         integer(c_size_t), value :: size
      end function cudaMallocHost

      integer(c_int) function cudaMallocManaged(devPtr, size, flags) &
         bind(C, name="cudaMallocManaged")
         import
         type(c_ptr), intent(out) :: devPtr
         integer(c_size_t), value :: size
         integer(c_int), value :: flags
      end function cudaMallocManaged

      integer(c_int) function cudaMallocMipmappedArray(mipmappedArray, desc, extent, numLevels, flags) &
         bind(C, name="cudaMallocMipmappedArray")
         import
         type(c_ptr), intent(out) :: mipmappedArray
         type(cudaChannelFormatDesc), intent(in) :: desc
         type(cudaExtent), value :: extent
         integer(c_int), value :: numLevels
         integer(c_int), value :: flags
      end function cudaMallocMipmappedArray

      integer(c_int) function cudaMallocPitch(devPtr, pitch, width, height) &
         bind(C, name="cudaMallocPitch")
         import
         type(c_ptr), intent(out) :: devPtr
         integer(c_size_t), intent(inout) :: pitch
         integer(c_size_t), value :: width
         integer(c_size_t), value :: height
      end function cudaMallocPitch

      integer(c_int) function cudaMemAdvise(devPtr, count, advice, device) &
         bind(C, name="cudaMemAdvise")
         import
         type(c_ptr), value :: devPtr
         integer(c_size_t), value :: count
         integer(c_int), value :: advice
         integer(c_int), value :: device
      end function cudaMemAdvise

      integer(c_int) function cudaMemAdvise_v2(devPtr, count, advice, location) &
         bind(C, name="cudaMemAdvise_v2")
         import
         type(c_ptr), value :: devPtr
         integer(c_size_t), value :: count
         integer(c_int), value :: advice
         type(cudaMemLocation), value :: location
      end function cudaMemAdvise_v2

      integer(c_int) function cudaMemGetInfo(free, total) &
         bind(C, name="cudaMemGetInfo")
         import
         integer(c_size_t), intent(inout) :: free
         integer(c_size_t), intent(inout) :: total
      end function cudaMemGetInfo

      integer(c_int) function cudaMemPoolCreate(memPool, poolProps) &
         bind(C, name="cudaMemPoolCreate")
         import
         type(c_ptr), intent(out) :: memPool
         type(cudaMemPoolProps), intent(in) :: poolProps
      end function cudaMemPoolCreate

      integer(c_int) function cudaMemPoolDestroy(memPool) &
         bind(C, name="cudaMemPoolDestroy")
         import
         type(c_ptr), value :: memPool
      end function cudaMemPoolDestroy

      integer(c_int) function cudaMemPoolExportPointer(exportData, ptr) &
         bind(C, name="cudaMemPoolExportPointer")
         import
         type(cudaMemPoolPtrExportData), intent(inout) :: exportData
         type(c_ptr), value :: ptr
      end function cudaMemPoolExportPointer

      integer(c_int) function cudaMemPoolExportToShareableHandle(shareableHandle, memPool, handleType, flags) &
         bind(C, name="cudaMemPoolExportToShareableHandle")
         import
         type(c_ptr), value :: shareableHandle
         type(c_ptr), value :: memPool
         integer(c_int), value :: handleType
         integer(c_int), value :: flags
      end function cudaMemPoolExportToShareableHandle

      integer(c_int) function cudaMemPoolGetAccess(flags, memPool, location) &
         bind(C, name="cudaMemPoolGetAccess")
         import
         integer(c_int), intent(out) :: flags
         type(c_ptr), value :: memPool
         type(cudaMemLocation), intent(inout) :: location
      end function cudaMemPoolGetAccess

      integer(c_int) function cudaMemPoolGetAttribute(memPool, attr, value) &
         bind(C, name="cudaMemPoolGetAttribute")
         import
         type(c_ptr), value :: memPool
         integer(c_int), value :: attr
         type(c_ptr), value :: value
      end function cudaMemPoolGetAttribute

      integer(c_int) function cudaMemPoolImportFromShareableHandle(memPool, shareableHandle, handleType, flags) &
         bind(C, name="cudaMemPoolImportFromShareableHandle")
         import
         type(c_ptr), intent(out) :: memPool
         type(c_ptr), value :: shareableHandle
         integer(c_int), value :: handleType
         integer(c_int), value :: flags
      end function cudaMemPoolImportFromShareableHandle

      integer(c_int) function cudaMemPoolImportPointer(ptr, memPool, exportData) &
         bind(C, name="cudaMemPoolImportPointer")
         import
         type(c_ptr), intent(out) :: ptr
         type(c_ptr), value :: memPool
         type(cudaMemPoolPtrExportData), intent(inout) :: exportData
      end function cudaMemPoolImportPointer

      integer(c_int) function cudaMemPoolSetAccess(memPool, descList, count) &
         bind(C, name="cudaMemPoolSetAccess")
         import
         type(c_ptr), value :: memPool
         type(cudaMemAccessDesc), intent(in) :: descList
         integer(c_size_t), value :: count
      end function cudaMemPoolSetAccess

      integer(c_int) function cudaMemPoolSetAttribute(memPool, attr, value) &
         bind(C, name="cudaMemPoolSetAttribute")
         import
         type(c_ptr), value :: memPool
         integer(c_int), value :: attr
         type(c_ptr), value :: value
      end function cudaMemPoolSetAttribute

      integer(c_int) function cudaMemPoolTrimTo(memPool, minBytesToKeep) &
         bind(C, name="cudaMemPoolTrimTo")
         import
         type(c_ptr), value :: memPool
         integer(c_size_t), value :: minBytesToKeep
      end function cudaMemPoolTrimTo

      integer(c_int) function cudaMemPrefetchAsync(devPtr, count, dstDevice, stream) &
         bind(C, name="cudaMemPrefetchAsync")
         import
         type(c_ptr), value :: devPtr
         integer(c_size_t), value :: count
         integer(c_int), value :: dstDevice
         type(c_ptr), value :: stream
      end function cudaMemPrefetchAsync

      integer(c_int) function cudaMemPrefetchAsync_v2(devPtr, count, location, flags, stream) &
         bind(C, name="cudaMemPrefetchAsync_v2")
         import
         type(c_ptr), value :: devPtr
         integer(c_size_t), value :: count
         type(cudaMemLocation), value :: location
         integer(c_int), value :: flags
         type(c_ptr), value :: stream
      end function cudaMemPrefetchAsync_v2

      integer(c_int) function cudaMemRangeGetAttribute(data, dataSize, attribute, devPtr, count) &
         bind(C, name="cudaMemRangeGetAttribute")
         import
         type(c_ptr), value :: data
         integer(c_size_t), value :: dataSize
         integer(c_int), value :: attribute
         type(c_ptr), value :: devPtr
         integer(c_size_t), value :: count
      end function cudaMemRangeGetAttribute

      integer(c_int) function cudaMemRangeGetAttributes(data, dataSizes, attributes, numAttributes, devPtr, count) &
         bind(C, name="cudaMemRangeGetAttributes")
         import
         type(c_ptr), intent(out) :: data
         integer(c_size_t), intent(inout) :: dataSizes
         integer(c_int), intent(out) :: attributes
         integer(c_size_t), value :: numAttributes
         type(c_ptr), value :: devPtr
         integer(c_size_t), value :: count
      end function cudaMemRangeGetAttributes

      integer(c_int) function cudaMemcpy(dst, src, count, kind) &
         bind(C, name="cudaMemcpy")
         import
         type(c_ptr), value :: dst
         type(c_ptr), value :: src
         integer(c_size_t), value :: count
         integer(c_int), value :: kind
      end function cudaMemcpy

      integer(c_int) function cudaMemcpy2D(dst, dpitch, src, spitch, width, height, kind) &
         bind(C, name="cudaMemcpy2D")
         import
         type(c_ptr), value :: dst
         integer(c_size_t), value :: dpitch
         type(c_ptr), value :: src
         integer(c_size_t), value :: spitch
         integer(c_size_t), value :: width
         integer(c_size_t), value :: height
         integer(c_int), value :: kind
      end function cudaMemcpy2D

      integer(c_int) function cudaMemcpy2DArrayToArray( &
         dst, wOffsetDst, hOffsetDst, src, wOffsetSrc, hOffsetSrc, width, height, kind) &
         bind(C, name="cudaMemcpy2DArrayToArray")
         import
         type(c_ptr), value :: dst
         integer(c_size_t), value :: wOffsetDst
         integer(c_size_t), value :: hOffsetDst
         type(c_ptr), value :: src
         integer(c_size_t), value :: wOffsetSrc
         integer(c_size_t), value :: hOffsetSrc
         integer(c_size_t), value :: width
         integer(c_size_t), value :: height
         integer(c_int), value :: kind
      end function cudaMemcpy2DArrayToArray

      integer(c_int) function cudaMemcpy2DAsync(dst, dpitch, src, spitch, width, height, kind, stream) &
         bind(C, name="cudaMemcpy2DAsync")
         import
         type(c_ptr), value :: dst
         integer(c_size_t), value :: dpitch
         type(c_ptr), value :: src
         integer(c_size_t), value :: spitch
         integer(c_size_t), value :: width
         integer(c_size_t), value :: height
         integer(c_int), value :: kind
         type(c_ptr), value :: stream
      end function cudaMemcpy2DAsync

      integer(c_int) function cudaMemcpy2DFromArray(dst, dpitch, src, wOffset, hOffset, width, height, kind) &
         bind(C, name="cudaMemcpy2DFromArray")
         import
         type(c_ptr), value :: dst
         integer(c_size_t), value :: dpitch
         type(c_ptr), value :: src
         integer(c_size_t), value :: wOffset
         integer(c_size_t), value :: hOffset
         integer(c_size_t), value :: width
         integer(c_size_t), value :: height
         integer(c_int), value :: kind
      end function cudaMemcpy2DFromArray

      integer(c_int) function cudaMemcpy2DFromArrayAsync( &
         dst, dpitch, src, wOffset, hOffset, width, height, kind, stream) &
         bind(C, name="cudaMemcpy2DFromArrayAsync")
         import
         type(c_ptr), value :: dst
         integer(c_size_t), value :: dpitch
         type(c_ptr), value :: src
         integer(c_size_t), value :: wOffset
         integer(c_size_t), value :: hOffset
         integer(c_size_t), value :: width
         integer(c_size_t), value :: height
         integer(c_int), value :: kind
         type(c_ptr), value :: stream
      end function cudaMemcpy2DFromArrayAsync

      integer(c_int) function cudaMemcpy2DToArray(dst, wOffset, hOffset, src, spitch, width, height, kind) &
         bind(C, name="cudaMemcpy2DToArray")
         import
         type(c_ptr), value :: dst
         integer(c_size_t), value :: wOffset
         integer(c_size_t), value :: hOffset
         type(c_ptr), value :: src
         integer(c_size_t), value :: spitch
         integer(c_size_t), value :: width
         integer(c_size_t), value :: height
         integer(c_int), value :: kind
      end function cudaMemcpy2DToArray

      integer(c_int) function cudaMemcpy2DToArrayAsync( &
         dst, wOffset, hOffset, src, spitch, width, height, kind, stream) &
         bind(C, name="cudaMemcpy2DToArrayAsync")
         import
         type(c_ptr), value :: dst
         integer(c_size_t), value :: wOffset
         integer(c_size_t), value :: hOffset
         type(c_ptr), value :: src
         integer(c_size_t), value :: spitch
         integer(c_size_t), value :: width
         integer(c_size_t), value :: height
         integer(c_int), value :: kind
         type(c_ptr), value :: stream
      end function cudaMemcpy2DToArrayAsync

      integer(c_int) function cudaMemcpy3D(p) &
         bind(C, name="cudaMemcpy3D")
         import
         type(cudaMemcpy3DParms), intent(in) :: p
      end function cudaMemcpy3D

      integer(c_int) function cudaMemcpy3DAsync(p, stream) &
         bind(C, name="cudaMemcpy3DAsync")
         import
         type(cudaMemcpy3DParms), intent(in) :: p
         type(c_ptr), value :: stream
      end function cudaMemcpy3DAsync

      integer(c_int) function cudaMemcpy3DBatchAsync(numOps, opList, failIdx, flags, stream) &
         bind(C, name="cudaMemcpy3DBatchAsync")
         import
         integer(c_size_t), value :: numOps
         type(cudaMemcpy3DBatchOp), intent(inout) :: opList
         integer(c_size_t), intent(inout) :: failIdx
         integer(c_long_long), value :: flags
         type(c_ptr), value :: stream
      end function cudaMemcpy3DBatchAsync

      integer(c_int) function cudaMemcpy3DPeer(p) &
         bind(C, name="cudaMemcpy3DPeer")
         import
         type(cudaMemcpy3DPeerParms), intent(in) :: p
      end function cudaMemcpy3DPeer

      integer(c_int) function cudaMemcpy3DPeerAsync(p, stream) &
         bind(C, name="cudaMemcpy3DPeerAsync")
         import
         type(cudaMemcpy3DPeerParms), intent(in) :: p
         type(c_ptr), value :: stream
      end function cudaMemcpy3DPeerAsync

      integer(c_int) function cudaMemcpyArrayToArray( &
         dst, wOffsetDst, hOffsetDst, src, wOffsetSrc, hOffsetSrc, count, kind) &
         bind(C, name="cudaMemcpyArrayToArray")
         import
         type(c_ptr), value :: dst
         integer(c_size_t), value :: wOffsetDst
         integer(c_size_t), value :: hOffsetDst
         type(c_ptr), value :: src
         integer(c_size_t), value :: wOffsetSrc
         integer(c_size_t), value :: hOffsetSrc
         integer(c_size_t), value :: count
         integer(c_int), value :: kind
      end function cudaMemcpyArrayToArray

      integer(c_int) function cudaMemcpyAsync(dst, src, count, kind, stream) &
         bind(C, name="cudaMemcpyAsync")
         import
         type(c_ptr), value :: dst
         type(c_ptr), value :: src
         integer(c_size_t), value :: count
         integer(c_int), value :: kind
         type(c_ptr), value :: stream
      end function cudaMemcpyAsync

      integer(c_int) function cudaMemcpyBatchAsync( &
         dsts, srcs, sizes, count, attrs, attrsIdxs, numAttrs, failIdx, stream) &
         bind(C, name="cudaMemcpyBatchAsync")
         import
         type(c_ptr), intent(out) :: dsts
         type(c_ptr), intent(out) :: srcs
         integer(c_size_t), intent(inout) :: sizes
         integer(c_size_t), value :: count
         type(cudaMemcpyAttributes), intent(inout) :: attrs
         integer(c_size_t), intent(inout) :: attrsIdxs
         integer(c_size_t), value :: numAttrs
         integer(c_size_t), intent(inout) :: failIdx
         type(c_ptr), value :: stream
      end function cudaMemcpyBatchAsync

      integer(c_int) function cudaMemcpyFromArray(dst, src, wOffset, hOffset, count, kind) &
         bind(C, name="cudaMemcpyFromArray")
         import
         type(c_ptr), value :: dst
         type(c_ptr), value :: src
         integer(c_size_t), value :: wOffset
         integer(c_size_t), value :: hOffset
         integer(c_size_t), value :: count
         integer(c_int), value :: kind
      end function cudaMemcpyFromArray

      integer(c_int) function cudaMemcpyFromArrayAsync(dst, src, wOffset, hOffset, count, kind, stream) &
         bind(C, name="cudaMemcpyFromArrayAsync")
         import
         type(c_ptr), value :: dst
         type(c_ptr), value :: src
         integer(c_size_t), value :: wOffset
         integer(c_size_t), value :: hOffset
         integer(c_size_t), value :: count
         integer(c_int), value :: kind
         type(c_ptr), value :: stream
      end function cudaMemcpyFromArrayAsync

      integer(c_int) function cudaMemcpyFromSymbol(dst, symbol, count, offset, kind) &
         bind(C, name="cudaMemcpyFromSymbol")
         import
         type(c_ptr), value :: dst
         type(c_ptr), value :: symbol
         integer(c_size_t), value :: count
         integer(c_size_t), value :: offset
         integer(c_int), value :: kind
      end function cudaMemcpyFromSymbol

      integer(c_int) function cudaMemcpyFromSymbolAsync(dst, symbol, count, offset, kind, stream) &
         bind(C, name="cudaMemcpyFromSymbolAsync")
         import
         type(c_ptr), value :: dst
         type(c_ptr), value :: symbol
         integer(c_size_t), value :: count
         integer(c_size_t), value :: offset
         integer(c_int), value :: kind
         type(c_ptr), value :: stream
      end function cudaMemcpyFromSymbolAsync

      integer(c_int) function cudaMemcpyPeer(dst, dstDevice, src, srcDevice, count) &
         bind(C, name="cudaMemcpyPeer")
         import
         type(c_ptr), value :: dst
         integer(c_int), value :: dstDevice
         type(c_ptr), value :: src
         integer(c_int), value :: srcDevice
         integer(c_size_t), value :: count
      end function cudaMemcpyPeer

      integer(c_int) function cudaMemcpyPeerAsync(dst, dstDevice, src, srcDevice, count, stream) &
         bind(C, name="cudaMemcpyPeerAsync")
         import
         type(c_ptr), value :: dst
         integer(c_int), value :: dstDevice
         type(c_ptr), value :: src
         integer(c_int), value :: srcDevice
         integer(c_size_t), value :: count
         type(c_ptr), value :: stream
      end function cudaMemcpyPeerAsync

      integer(c_int) function cudaMemcpyToArray(dst, wOffset, hOffset, src, count, kind) &
         bind(C, name="cudaMemcpyToArray")
         import
         type(c_ptr), value :: dst
         integer(c_size_t), value :: wOffset
         integer(c_size_t), value :: hOffset
         type(c_ptr), value :: src
         integer(c_size_t), value :: count
         integer(c_int), value :: kind
      end function cudaMemcpyToArray

      integer(c_int) function cudaMemcpyToArrayAsync(dst, wOffset, hOffset, src, count, kind, stream) &
         bind(C, name="cudaMemcpyToArrayAsync")
         import
         type(c_ptr), value :: dst
         integer(c_size_t), value :: wOffset
         integer(c_size_t), value :: hOffset
         type(c_ptr), value :: src
         integer(c_size_t), value :: count
         integer(c_int), value :: kind
         type(c_ptr), value :: stream
      end function cudaMemcpyToArrayAsync

      integer(c_int) function cudaMemcpyToSymbol(symbol, src, count, offset, kind) &
         bind(C, name="cudaMemcpyToSymbol")
         import
         type(c_ptr), value :: symbol
         type(c_ptr), value :: src
         integer(c_size_t), value :: count
         integer(c_size_t), value :: offset
         integer(c_int), value :: kind
      end function cudaMemcpyToSymbol

      integer(c_int) function cudaMemcpyToSymbolAsync(symbol, src, count, offset, kind, stream) &
         bind(C, name="cudaMemcpyToSymbolAsync")
         import
         type(c_ptr), value :: symbol
         type(c_ptr), value :: src
         integer(c_size_t), value :: count
         integer(c_size_t), value :: offset
         integer(c_int), value :: kind
         type(c_ptr), value :: stream
      end function cudaMemcpyToSymbolAsync

      integer(c_int) function cudaMemset(devPtr, value, count) &
         bind(C, name="cudaMemset")
         import
         type(c_ptr), value :: devPtr
         integer(c_int), value :: value
         integer(c_size_t), value :: count
      end function cudaMemset

      integer(c_int) function cudaMemset2D(devPtr, pitch, value, width, height) &
         bind(C, name="cudaMemset2D")
         import
         type(c_ptr), value :: devPtr
         integer(c_size_t), value :: pitch
         integer(c_int), value :: value
         integer(c_size_t), value :: width
         integer(c_size_t), value :: height
      end function cudaMemset2D

      integer(c_int) function cudaMemset2DAsync(devPtr, pitch, value, width, height, stream) &
         bind(C, name="cudaMemset2DAsync")
         import
         type(c_ptr), value :: devPtr
         integer(c_size_t), value :: pitch
         integer(c_int), value :: value
         integer(c_size_t), value :: width
         integer(c_size_t), value :: height
         type(c_ptr), value :: stream
      end function cudaMemset2DAsync

      integer(c_int) function cudaMemset3D(pitchedDevPtr, value, extent) &
         bind(C, name="cudaMemset3D")
         import
         type(cudaPitchedPtr), value :: pitchedDevPtr
         integer(c_int), value :: value
         type(cudaExtent), value :: extent
      end function cudaMemset3D

      integer(c_int) function cudaMemset3DAsync(pitchedDevPtr, value, extent, stream) &
         bind(C, name="cudaMemset3DAsync")
         import
         type(cudaPitchedPtr), value :: pitchedDevPtr
         integer(c_int), value :: value
         type(cudaExtent), value :: extent
         type(c_ptr), value :: stream
      end function cudaMemset3DAsync

      integer(c_int) function cudaMemsetAsync(devPtr, value, count, stream) &
         bind(C, name="cudaMemsetAsync")
         import
         type(c_ptr), value :: devPtr
         integer(c_int), value :: value
         integer(c_size_t), value :: count
         type(c_ptr), value :: stream
      end function cudaMemsetAsync

      integer(c_int) function cudaMipmappedArrayGetMemoryRequirements(memoryRequirements, mipmap, device) &
         bind(C, name="cudaMipmappedArrayGetMemoryRequirements")
         import
         type(cudaArrayMemoryRequirements), intent(inout) :: memoryRequirements
         type(c_ptr), value :: mipmap
         integer(c_int), value :: device
      end function cudaMipmappedArrayGetMemoryRequirements

      integer(c_int) function cudaMipmappedArrayGetSparseProperties(sparseProperties, mipmap) &
         bind(C, name="cudaMipmappedArrayGetSparseProperties")
         import
         type(cudaArraySparseProperties), intent(inout) :: sparseProperties
         type(c_ptr), value :: mipmap
      end function cudaMipmappedArrayGetSparseProperties

      integer(c_int) function cudaOccupancyAvailableDynamicSMemPerBlock(dynamicSmemSize, func, numBlocks, blockSize) &
         bind(C, name="cudaOccupancyAvailableDynamicSMemPerBlock")
         import
         integer(c_size_t), intent(inout) :: dynamicSmemSize
         type(c_ptr), value :: func
         integer(c_int), value :: numBlocks
         integer(c_int), value :: blockSize
      end function cudaOccupancyAvailableDynamicSMemPerBlock

      integer(c_int) function cudaOccupancyMaxActiveBlocksPerMultiprocessor( &
         numBlocks, func, blockSize, dynamicSMemSize) &
         bind(C, name="cudaOccupancyMaxActiveBlocksPerMultiprocessor")
         import
         integer(c_int), intent(inout) :: numBlocks
         type(c_ptr), value :: func
         integer(c_int), value :: blockSize
         integer(c_size_t), value :: dynamicSMemSize
      end function cudaOccupancyMaxActiveBlocksPerMultiprocessor

      integer(c_int) function cudaOccupancyMaxActiveBlocksPerMultiprocessorWithFlags( &
         numBlocks, func, blockSize, dynamicSMemSize, flags) &
         bind(C, name="cudaOccupancyMaxActiveBlocksPerMultiprocessorWithFlags")
         import
         integer(c_int), intent(inout) :: numBlocks
         type(c_ptr), value :: func
         integer(c_int), value :: blockSize
         integer(c_size_t), value :: dynamicSMemSize
         integer(c_int), value :: flags
      end function cudaOccupancyMaxActiveBlocksPerMultiprocessorWithFlags

      integer(c_int) function cudaOccupancyMaxActiveClusters(numClusters, func, launchConfig) &
         bind(C, name="cudaOccupancyMaxActiveClusters")
         import
         integer(c_int), intent(inout) :: numClusters
         type(c_ptr), value :: func
         type(cudaLaunchConfig_st), intent(in) :: launchConfig
      end function cudaOccupancyMaxActiveClusters

      integer(c_int) function cudaOccupancyMaxPotentialClusterSize(clusterSize, func, launchConfig) &
         bind(C, name="cudaOccupancyMaxPotentialClusterSize")
         import
         integer(c_int), intent(inout) :: clusterSize
         type(c_ptr), value :: func
         type(cudaLaunchConfig_st), intent(in) :: launchConfig
      end function cudaOccupancyMaxPotentialClusterSize

      integer(c_int) function cudaPeekAtLastError() &
         bind(C, name="cudaPeekAtLastError")
         import
      end function cudaPeekAtLastError

      integer(c_int) function cudaPointerGetAttributes(attributes, ptr) &
         bind(C, name="cudaPointerGetAttributes")
         import
         type(cudaPointerAttributes), intent(inout) :: attributes
         type(c_ptr), value :: ptr
      end function cudaPointerGetAttributes

      integer(c_int) function cudaProfilerStart() &
         bind(C, name="cudaProfilerStart")
         import
      end function cudaProfilerStart

      integer(c_int) function cudaProfilerStop() &
         bind(C, name="cudaProfilerStop")
         import
      end function cudaProfilerStop

      integer(c_int) function cudaRuntimeGetVersion(runtimeVersion) &
         bind(C, name="cudaRuntimeGetVersion")
         import
         integer(c_int), intent(inout) :: runtimeVersion
      end function cudaRuntimeGetVersion

      integer(c_int) function cudaSetDevice(device) &
         bind(C, name="cudaSetDevice")
         import
         integer(c_int), value :: device
      end function cudaSetDevice

      integer(c_int) function cudaSetDeviceFlags(flags) &
         bind(C, name="cudaSetDeviceFlags")
         import
         integer(c_int), value :: flags
      end function cudaSetDeviceFlags

      integer(c_int) function cudaSetDoubleForDevice(d) &
         bind(C, name="cudaSetDoubleForDevice")
         import
         real(c_double), intent(inout) :: d
      end function cudaSetDoubleForDevice

      integer(c_int) function cudaSetDoubleForHost(d) &
         bind(C, name="cudaSetDoubleForHost")
         import
         real(c_double), intent(inout) :: d
      end function cudaSetDoubleForHost

      integer(c_int) function cudaSetValidDevices(device_arr, len) &
         bind(C, name="cudaSetValidDevices")
         import
         integer(c_int), intent(inout) :: device_arr
         integer(c_int), value :: len
      end function cudaSetValidDevices

      integer(c_int) function cudaSignalExternalSemaphoresAsync_v2(extSemArray, paramsArray, numExtSems, stream) &
         bind(C, name="cudaSignalExternalSemaphoresAsync_v2")
         import
         type(c_ptr), intent(out) :: extSemArray
         type(cudaExternalSemaphoreSignalParams), intent(in) :: paramsArray
         integer(c_int), value :: numExtSems
         type(c_ptr), value :: stream
      end function cudaSignalExternalSemaphoresAsync_v2

      integer(c_int) function cudaStreamAddCallback(stream, callback, userData, flags) &
         bind(C, name="cudaStreamAddCallback")
         import
         type(c_ptr), value :: stream
         type(c_funptr), value :: callback
         type(c_ptr), value :: userData
         integer(c_int), value :: flags
      end function cudaStreamAddCallback

      integer(c_int) function cudaStreamAttachMemAsync(stream, devPtr, length, flags) &
         bind(C, name="cudaStreamAttachMemAsync")
         import
         type(c_ptr), value :: stream
         type(c_ptr), value :: devPtr
         integer(c_size_t), value :: length
         integer(c_int), value :: flags
      end function cudaStreamAttachMemAsync

      integer(c_int) function cudaStreamBeginCapture(stream, mode) &
         bind(C, name="cudaStreamBeginCapture")
         import
         type(c_ptr), value :: stream
         integer(c_int), value :: mode
      end function cudaStreamBeginCapture

      integer(c_int) function cudaStreamBeginCaptureToGraph( &
         stream, graph, dependencies, dependencyData, numDependencies, mode) &
         bind(C, name="cudaStreamBeginCaptureToGraph")
         import
         type(c_ptr), value :: stream
         type(c_ptr), value :: graph
         type(c_ptr), intent(out) :: dependencies
         type(cudaGraphEdgeData_st), intent(in) :: dependencyData
         integer(c_size_t), value :: numDependencies
         integer(c_int), value :: mode
      end function cudaStreamBeginCaptureToGraph

      integer(c_int) function cudaStreamCopyAttributes(dst, src) &
         bind(C, name="cudaStreamCopyAttributes")
         import
         type(c_ptr), value :: dst
         type(c_ptr), value :: src
      end function cudaStreamCopyAttributes

      integer(c_int) function cudaStreamCreate(pStream) &
         bind(C, name="cudaStreamCreate")
         import
         type(c_ptr), intent(out) :: pStream
      end function cudaStreamCreate

      integer(c_int) function cudaStreamCreateWithFlags(pStream, flags) &
         bind(C, name="cudaStreamCreateWithFlags")
         import
         type(c_ptr), intent(out) :: pStream
         integer(c_int), value :: flags
      end function cudaStreamCreateWithFlags

      integer(c_int) function cudaStreamCreateWithPriority(pStream, flags, priority) &
         bind(C, name="cudaStreamCreateWithPriority")
         import
         type(c_ptr), intent(out) :: pStream
         integer(c_int), value :: flags
         integer(c_int), value :: priority
      end function cudaStreamCreateWithPriority

      integer(c_int) function cudaStreamDestroy(stream) &
         bind(C, name="cudaStreamDestroy")
         import
         type(c_ptr), value :: stream
      end function cudaStreamDestroy

      integer(c_int) function cudaStreamEndCapture(stream, pGraph) &
         bind(C, name="cudaStreamEndCapture")
         import
         type(c_ptr), value :: stream
         type(c_ptr), intent(out) :: pGraph
      end function cudaStreamEndCapture

      integer(c_int) function cudaStreamGetAttribute(hStream, attr, value_out) &
         bind(C, name="cudaStreamGetAttribute")
         import
         type(c_ptr), value :: hStream
         integer(c_int), value :: attr
         type(cudaLaunchAttributeValue), intent(inout) :: value_out
      end function cudaStreamGetAttribute

      integer(c_int) function cudaStreamGetCaptureInfo_v2( &
         stream, captureStatus_out, id_out, graph_out, dependencies_out, numDependencies_out) &
         bind(C, name="cudaStreamGetCaptureInfo_v2")
         import
         type(c_ptr), value :: stream
         integer(c_int), intent(out) :: captureStatus_out
         integer(c_long_long), intent(inout) :: id_out
         type(c_ptr), intent(out) :: graph_out
         type(c_ptr), intent(out) :: dependencies_out
         integer(c_size_t), intent(inout) :: numDependencies_out
      end function cudaStreamGetCaptureInfo_v2

      integer(c_int) function cudaStreamGetCaptureInfo_v3( &
         stream, captureStatus_out, id_out, graph_out, dependencies_out, edgeData_out, numDependencies_out) &
         bind(C, name="cudaStreamGetCaptureInfo_v3")
         import
         type(c_ptr), value :: stream
         integer(c_int), intent(out) :: captureStatus_out
         integer(c_long_long), intent(inout) :: id_out
         type(c_ptr), intent(out) :: graph_out
         type(c_ptr), intent(out) :: dependencies_out
         type(cudaGraphEdgeData_st), intent(in) :: edgeData_out
         integer(c_size_t), intent(inout) :: numDependencies_out
      end function cudaStreamGetCaptureInfo_v3

      integer(c_int) function cudaStreamGetDevice(hStream, device) &
         bind(C, name="cudaStreamGetDevice")
         import
         type(c_ptr), value :: hStream
         integer(c_int), intent(inout) :: device
      end function cudaStreamGetDevice

      integer(c_int) function cudaStreamGetFlags(hStream, flags) &
         bind(C, name="cudaStreamGetFlags")
         import
         type(c_ptr), value :: hStream
         integer(c_int), intent(inout) :: flags
      end function cudaStreamGetFlags

      integer(c_int) function cudaStreamGetId(hStream, streamId) &
         bind(C, name="cudaStreamGetId")
         import
         type(c_ptr), value :: hStream
         integer(c_long_long), intent(inout) :: streamId
      end function cudaStreamGetId

      integer(c_int) function cudaStreamGetPriority(hStream, priority) &
         bind(C, name="cudaStreamGetPriority")
         import
         type(c_ptr), value :: hStream
         integer(c_int), intent(inout) :: priority
      end function cudaStreamGetPriority

      integer(c_int) function cudaStreamIsCapturing(stream, pCaptureStatus) &
         bind(C, name="cudaStreamIsCapturing")
         import
         type(c_ptr), value :: stream
         integer(c_int), intent(out) :: pCaptureStatus
      end function cudaStreamIsCapturing

      integer(c_int) function cudaStreamQuery(stream) &
         bind(C, name="cudaStreamQuery")
         import
         type(c_ptr), value :: stream
      end function cudaStreamQuery

      integer(c_int) function cudaStreamSetAttribute(hStream, attr, value) &
         bind(C, name="cudaStreamSetAttribute")
         import
         type(c_ptr), value :: hStream
         integer(c_int), value :: attr
         type(cudaLaunchAttributeValue), intent(in) :: value
      end function cudaStreamSetAttribute

      integer(c_int) function cudaStreamSynchronize(stream) &
         bind(C, name="cudaStreamSynchronize")
         import
         type(c_ptr), value :: stream
      end function cudaStreamSynchronize

      integer(c_int) function cudaStreamUpdateCaptureDependencies(stream, dependencies, numDependencies, flags) &
         bind(C, name="cudaStreamUpdateCaptureDependencies")
         import
         type(c_ptr), value :: stream
         type(c_ptr), intent(out) :: dependencies
         integer(c_size_t), value :: numDependencies
         integer(c_int), value :: flags
      end function cudaStreamUpdateCaptureDependencies

      integer(c_int) function cudaStreamUpdateCaptureDependencies_v2( &
         stream, dependencies, dependencyData, numDependencies, flags) &
         bind(C, name="cudaStreamUpdateCaptureDependencies_v2")
         import
         type(c_ptr), value :: stream
         type(c_ptr), intent(out) :: dependencies
         type(cudaGraphEdgeData_st), intent(in) :: dependencyData
         integer(c_size_t), value :: numDependencies
         integer(c_int), value :: flags
      end function cudaStreamUpdateCaptureDependencies_v2

      integer(c_int) function cudaStreamWaitEvent(stream, event, flags) &
         bind(C, name="cudaStreamWaitEvent")
         import
         type(c_ptr), value :: stream
         type(c_ptr), value :: event
         integer(c_int), value :: flags
      end function cudaStreamWaitEvent

      integer(c_int) function cudaThreadExchangeStreamCaptureMode(mode) &
         bind(C, name="cudaThreadExchangeStreamCaptureMode")
         import
         integer(c_int), intent(out) :: mode
      end function cudaThreadExchangeStreamCaptureMode

      integer(c_int) function cudaThreadExit() &
         bind(C, name="cudaThreadExit")
         import
      end function cudaThreadExit

      integer(c_int) function cudaThreadGetCacheConfig(pCacheConfig) &
         bind(C, name="cudaThreadGetCacheConfig")
         import
         integer(c_int), intent(out) :: pCacheConfig
      end function cudaThreadGetCacheConfig

      integer(c_int) function cudaThreadGetLimit(pValue, limit) &
         bind(C, name="cudaThreadGetLimit")
         import
         integer(c_size_t), intent(inout) :: pValue
         integer(c_int), value :: limit
      end function cudaThreadGetLimit

      integer(c_int) function cudaThreadSetCacheConfig(cacheConfig) &
         bind(C, name="cudaThreadSetCacheConfig")
         import
         integer(c_int), value :: cacheConfig
      end function cudaThreadSetCacheConfig

      integer(c_int) function cudaThreadSetLimit(limit, value) &
         bind(C, name="cudaThreadSetLimit")
         import
         integer(c_int), value :: limit
         integer(c_size_t), value :: value
      end function cudaThreadSetLimit

      integer(c_int) function cudaThreadSynchronize() &
         bind(C, name="cudaThreadSynchronize")
         import
      end function cudaThreadSynchronize

      integer(c_int) function cudaUserObjectCreate(object_out, ptr, destroy, initialRefcount, flags) &
         bind(C, name="cudaUserObjectCreate")
         import
         type(c_ptr), intent(out) :: object_out
         type(c_ptr), value :: ptr
         type(c_funptr), value :: destroy
         integer(c_int), value :: initialRefcount
         integer(c_int), value :: flags
      end function cudaUserObjectCreate

      integer(c_int) function cudaUserObjectRelease(object, count) &
         bind(C, name="cudaUserObjectRelease")
         import
         type(c_ptr), value :: object
         integer(c_int), value :: count
      end function cudaUserObjectRelease

      integer(c_int) function cudaUserObjectRetain(object, count) &
         bind(C, name="cudaUserObjectRetain")
         import
         type(c_ptr), value :: object
         integer(c_int), value :: count
      end function cudaUserObjectRetain

      integer(c_int) function cudaWaitExternalSemaphoresAsync_v2(extSemArray, paramsArray, numExtSems, stream) &
         bind(C, name="cudaWaitExternalSemaphoresAsync_v2")
         import
         type(c_ptr), intent(out) :: extSemArray
         type(cudaExternalSemaphoreWaitParams), intent(in) :: paramsArray
         integer(c_int), value :: numExtSems
         type(c_ptr), value :: stream
      end function cudaWaitExternalSemaphoresAsync_v2

   end interface

end module cuda_runtime
