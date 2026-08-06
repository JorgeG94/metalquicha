!! Per-rank cuEST handle and GPU device binding
module mqc_cuest_context
   !! Owns the cuEST library handle for this MPI rank.
   !!
   !! Handle creation is expensive and the handle is reusable across any number
   !! of molecules, so it is created once per rank and shared by every fragment
   !! that rank evaluates. This matters for fragmented runs, where a rank may
   !! process thousands of fragments.
   !!
   !! GPU selection is by *node-local* rank: with several ranks per node the
   !! devices are handed out round-robin, so ranks 0,1,2,3 on a 4-GPU node land
   !! on distinct devices. Oversubscription (more ranks than GPUs) is allowed
   !! and simply shares devices.
   !!
   !! cuEST ships GPU code for sm_80 and newer only. On anything older --
   !! including Volta/V100 -- `cuestCreate` returns
   !! CUEST_STATUS_UNSUPPORTED_ARCHITECTURE, which surfaces here as an error
   !! rather than a crash.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_int, c_int64_t, c_associated
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_cuest_runtime, only: cuest_status_check, cuda_status_check, &
                                cublas_status_check, device_alloc, device_free
   use cuest, only: cuestCreate, cuestDestroy, cuestParametersCreate, &
                    cuestParametersDestroy, cuestSetMathMode, &
                    CUEST_HANDLE_PARAMETERS, CUEST_NATIVE_FP64_MATH_MODE
   use cublas, only: cublasCreate, cublasDestroy
   use cuda_runtime, only: cudaSetDevice, cudaGetDeviceCount
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   implicit none
   private

   public :: cuest_context_t        !! Per-rank cuEST handle
   public :: get_cuest_context      !! Lazily initialized process-wide context
   public :: finalize_cuest_context  !! Release the process-wide context

   type :: device_pool_t
      !! A device buffer that grows to a high-water mark and is never shrunk
      !!
      !! Fragmented runs evaluate thousands of small molecules on one rank. A
      !! `cudaMalloc` per fragment would dominate the runtime -- it is a
      !! synchronizing call -- so scratch is allocated once at the largest size
      !! any fragment has needed and handed out again for every fragment that
      !! fits. Growth reallocates; there is no shrink, by design.
      type(c_ptr) :: ptr = c_null_ptr        !! Device address, null until first use
      integer(c_int64_t) :: capacity = 0     !! Current capacity in doubles
   contains
      procedure :: ensure => pool_ensure     !! Grow to hold at least n doubles
      procedure :: release => pool_release   !! Free the buffer
   end type device_pool_t

   type :: cuest_context_t
      !! A live cuEST handle bound to one GPU, plus reusable device scratch
      type(c_ptr) :: handle = c_null_ptr  !! cuestHandle_t
      type(c_ptr) :: cublas_handle = c_null_ptr
         !! cublasHandle_t, for the SCF's own linear algebra on device-resident
         !! matrices. Created alongside the cuEST handle because it is bound to
         !! the same device and has the same per-rank lifetime.
      integer :: device_id = -1           !! CUDA device this handle is bound to
      logical :: initialized = .false.    !! Whether `handle` is live

      ! Scratch shared by every fragment this rank evaluates. Fragment-local
      ! cuEST objects (basis, pair list, plans) are inherently geometry- and
      ! basis-specific and must still be rebuilt per fragment; these plain
      ! matrix buffers need not be.
      type(device_pool_t) :: scratch_density  !! Density / MO coefficient input
      type(device_pool_t) :: scratch_result   !! Matrix being built
      type(device_pool_t) :: scratch_c_occ    !! Occupied MO coefficients (alpha)
      type(device_pool_t) :: scratch_c_occ_beta  !! Beta occupied MOs (UKS)
      type(device_pool_t) :: scratch_result_beta  !! Second output matrix (UKS)
      type(device_pool_t) :: scratch_xyz      !! Atom coordinates
      type(device_pool_t) :: scratch_charges  !! Nuclear charges
      type(device_pool_t) :: scratch_gradient  !! natom x 3 gradient output
      type(device_pool_t) :: scratch_charge_gradient  !! Hellmann-Feynman half

      ! One output buffer per Fock term. `scratch_result` alone cannot serve
      ! all three: J, K and Vxc would each have to be pulled to the host before
      ! the next call overwrote it, which is what forced the round trip the
      ! device-resident path exists to remove.
      type(device_pool_t) :: scratch_j     !! Coulomb output
      type(device_pool_t) :: scratch_k     !! Exchange output
      type(device_pool_t) :: scratch_xc    !! XC potential output
      type(device_pool_t) :: scratch_fock  !! Assembled Fock, stays resident
      type(device_pool_t) :: scratch_core  !! Core Hamiltonian, uploaded once
      type(device_pool_t) :: scratch_ovlp  !! Overlap, uploaded once
   contains
      procedure :: create => context_create    !! Bind a device and create the handle
      procedure :: destroy => context_destroy  !! Release the handle and scratch
   end type cuest_context_t

   type(cuest_context_t), save, target :: shared_context
      !! Process-wide context, created on first use by `get_cuest_context`.
      !! One per process means one per MPI rank: each rank owns its handle,
      !! its device binding and its scratch pools, with nothing shared between
      !! ranks and so nothing to synchronize.

contains

   subroutine pool_ensure(this, n_doubles, label, error)
      !! Make sure the pool can hold `n_doubles` doubles, growing if it cannot
      class(device_pool_t), intent(inout) :: this
      integer(c_int64_t), intent(in) :: n_doubles
      character(len=*), intent(in) :: label
      type(error_t), intent(inout) :: error

      if (n_doubles <= this%capacity) return

      call device_free(this%ptr)
      this%capacity = 0
      call device_alloc(this%ptr, n_doubles, label, error)
      if (error%has_error()) return
      this%capacity = n_doubles
   end subroutine pool_ensure

   subroutine pool_release(this)
      !! Free the pooled buffer
      class(device_pool_t), intent(inout) :: this

      call device_free(this%ptr)
      this%capacity = 0
   end subroutine pool_release

   subroutine context_create(this, local_rank, error)
      !! Bind this rank to a GPU and create its cuEST handle
      class(cuest_context_t), intent(inout) :: this
      integer, intent(in) :: local_rank  !! Node-local MPI rank, 0-based
      type(error_t), intent(inout) :: error

      type(c_ptr) :: handle_params
      integer(c_int) :: device_count, status

      if (this%initialized) return

      ! ---- pick a device ---------------------------------------------------
      device_count = 0
      call cuda_status_check(cudaGetDeviceCount(device_count), &
                             "cudaGetDeviceCount", error)
      if (error%has_error()) return

      if (device_count <= 0) then
         call error%set(ERROR_VALIDATION, "No CUDA devices visible to this rank")
         return
      end if

      this%device_id = mod(local_rank, int(device_count))
      call cuda_status_check(cudaSetDevice(int(this%device_id, c_int)), &
                             "cudaSetDevice", error)
      if (error%has_error()) return

      ! ---- create the handle ----------------------------------------------
      handle_params = c_null_ptr
      call cuest_status_check(cuestParametersCreate(CUEST_HANDLE_PARAMETERS, handle_params), &
                              "cuestParametersCreate(handle)", error)
      if (error%has_error()) return

      call cuest_status_check(cuestCreate(handle_params, this%handle), &
                              "cuestCreate", error)

      status = cuestParametersDestroy(CUEST_HANDLE_PARAMETERS, handle_params)
      call cuest_status_check(status, "cuestParametersDestroy(handle)", error)
      if (error%has_error()) then
         call error%add_context("mqc_cuest_context:create (needs a GPU of sm_80 or newer)")
         return
      end if

      ! Fragment energies are differenced against each other in the MBE, so the
      ! integrals must be reproducible to full double precision -- no TF32.
      call cuest_status_check(cuestSetMathMode(this%handle, CUEST_NATIVE_FP64_MATH_MODE), &
                              "cuestSetMathMode(native FP64)", error)
      if (error%has_error()) return

      ! cuBLAS binds to whatever device is current, which cudaSetDevice above
      ! has already fixed. The handle is left in its default
      ! CUBLAS_POINTER_MODE_HOST, which is what the bindings assume.
      call cublas_status_check(cublasCreate(this%cublas_handle), "cublasCreate", error)
      if (error%has_error()) return

      this%initialized = .true.

      ! Say which device this rank took. Without it a misbinding -- every rank
      ! piling onto GPU 0 -- looks exactly like a slow run, and that is an
      ! expensive thing to diagnose after the fact.
      call logger%info("cuEST: node-local rank "//to_char(local_rank)// &
                       " bound to GPU "//to_char(this%device_id)// &
                       " of "//to_char(int(device_count)))
   end subroutine context_create

   subroutine context_destroy(this)
      !! Release the pooled scratch and the cuEST handle
      class(cuest_context_t), intent(inout) :: this
      integer(c_int) :: status

      call this%scratch_density%release()
      call this%scratch_result%release()
      call this%scratch_c_occ%release()
      call this%scratch_c_occ_beta%release()
      call this%scratch_result_beta%release()
      call this%scratch_xyz%release()
      call this%scratch_charges%release()
      call this%scratch_gradient%release()
      call this%scratch_charge_gradient%release()
      call this%scratch_j%release()
      call this%scratch_k%release()
      call this%scratch_xc%release()
      call this%scratch_fock%release()
      call this%scratch_core%release()
      call this%scratch_ovlp%release()

      if (c_associated(this%cublas_handle)) status = cublasDestroy(this%cublas_handle)
      this%cublas_handle = c_null_ptr

      if (c_associated(this%handle)) status = cuestDestroy(this%handle)
      this%handle = c_null_ptr
      this%device_id = -1
      this%initialized = .false.
   end subroutine context_destroy

   subroutine get_cuest_context(local_rank, context, error)
      !! Return the process-wide context, creating it on first call
      integer, intent(in) :: local_rank  !! Node-local MPI rank, 0-based
      type(cuest_context_t), pointer, intent(out) :: context
      type(error_t), intent(inout) :: error

      if (.not. shared_context%initialized) then
         call shared_context%create(local_rank, error)
         if (error%has_error()) then
            context => null()
            return
         end if
      end if
      context => shared_context
   end subroutine get_cuest_context

   subroutine finalize_cuest_context()
      !! Release the process-wide context, if one was ever created
      call shared_context%destroy()
   end subroutine finalize_cuest_context

end module mqc_cuest_context
