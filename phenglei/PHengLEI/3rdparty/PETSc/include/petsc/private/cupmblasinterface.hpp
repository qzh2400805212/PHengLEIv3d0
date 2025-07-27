#ifndef PETSCCUPMBLASINTERFACE_HPP
#define PETSCCUPMBLASINTERFACE_HPP

#if defined(__cplusplus)
  #include <petsc/private/cupminterface.hpp>
  #include <petsc/private/petscadvancedmacros.h>

namespace Petsc
{

namespace device
{

namespace cupm
{

namespace impl
{

  #define PetscCallCUPMBLAS(...) \
    do { \
      const cupmBlasError_t cberr_p_ = __VA_ARGS__; \
      if (PetscUnlikely(cberr_p_ != CUPMBLAS_STATUS_SUCCESS)) { \
        if (((cberr_p_ == CUPMBLAS_STATUS_NOT_INITIALIZED) || (cberr_p_ == CUPMBLAS_STATUS_ALLOC_FAILED)) && PetscDeviceInitialized(PETSC_DEVICE_CUPM())) { \
          SETERRQ(PETSC_COMM_SELF, PETSC_ERR_GPU_RESOURCE, \
                  "%s error %d (%s). Reports not initialized or alloc failed; " \
                  "this indicates the GPU may have run out resources", \
                  cupmBlasName(), static_cast<PetscErrorCode>(cberr_p_), cupmBlasGetErrorName(cberr_p_)); \
        } \
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_GPU, "%s error %d (%s)", cupmBlasName(), static_cast<PetscErrorCode>(cberr_p_), cupmBlasGetErrorName(cberr_p_)); \
      } \
    } while (0)

  // given cupmBlas<T>axpy() then
  // T = PETSC_CUPBLAS_FP_TYPE
  // given cupmBlas<T><u>nrm2() then
  // T = PETSC_CUPMBLAS_FP_INPUT_TYPE
  // u = PETSC_CUPMBLAS_FP_RETURN_TYPE
  #if PetscDefined(USE_COMPLEX)
    #if PetscDefined(USE_REAL_SINGLE)
      #define PETSC_CUPMBLAS_FP_TYPE_U       C
      #define PETSC_CUPMBLAS_FP_TYPE_L       c
      #define PETSC_CUPMBLAS_FP_INPUT_TYPE_U S
      #define PETSC_CUPMBLAS_FP_INPUT_TYPE_L s
    #elif PetscDefined(USE_REAL_DOUBLE)
      #define PETSC_CUPMBLAS_FP_TYPE_U       Z
      #define PETSC_CUPMBLAS_FP_TYPE_L       z
      #define PETSC_CUPMBLAS_FP_INPUT_TYPE_U D
      #define PETSC_CUPMBLAS_FP_INPUT_TYPE_L d
    #endif
    #define PETSC_CUPMBLAS_FP_RETURN_TYPE_U PETSC_CUPMBLAS_FP_TYPE_U
    #define PETSC_CUPMBLAS_FP_RETURN_TYPE_L PETSC_CUPMBLAS_FP_TYPE_L
  #else
    #if PetscDefined(USE_REAL_SINGLE)
      #define PETSC_CUPMBLAS_FP_TYPE_U S
      #define PETSC_CUPMBLAS_FP_TYPE_L s
    #elif PetscDefined(USE_REAL_DOUBLE)
      #define PETSC_CUPMBLAS_FP_TYPE_U D
      #define PETSC_CUPMBLAS_FP_TYPE_L d
    #endif
    #define PETSC_CUPMBLAS_FP_INPUT_TYPE_U PETSC_CUPMBLAS_FP_TYPE_U
    #define PETSC_CUPMBLAS_FP_INPUT_TYPE_L PETSC_CUPMBLAS_FP_TYPE_L
    #define PETSC_CUPMBLAS_FP_RETURN_TYPE_U
    #define PETSC_CUPMBLAS_FP_RETURN_TYPE_L
  #endif // USE_COMPLEX

  #if !defined(PETSC_CUPMBLAS_FP_TYPE_U) && !PetscDefined(USE_REAL___FLOAT128)
    #error "Unsupported floating-point type for CUDA/HIP BLAS"
  #endif

  // PETSC_CUPMBLAS_BUILD_BLAS_FUNCTION_ALIAS_MODIFIED() - Helper macro to build a "modified"
  // blas function whose return type does not match the input type
  //
  // input param:
  // func - base suffix of the blas function, e.g. nrm2
  //
  // notes:
  // requires PETSC_CUPMBLAS_FP_INPUT_TYPE to be defined as the blas floating point input type
  // letter ("S" for real/complex single, "D" for real/complex double).
  //
  // requires PETSC_CUPMBLAS_FP_RETURN_TYPE to be defined as the blas floating point output type
  // letter ("c" for complex single, "z" for complex double and <absolutely nothing> for real
  // single/double).
  //
  // In their infinite wisdom nvidia/amd have made the upper-case vs lower-case scheme
  // infuriatingly inconsistent...
  //
  // example usage:
  // #define PETSC_CUPMBLAS_FP_INPUT_TYPE  S
  // #define PETSC_CUPMBLAS_FP_RETURN_TYPE
  // PETSC_CUPMBLAS_BUILD_BLAS_FUNCTION_ALIAS_MODIFIED(nrm2) -> Snrm2
  //
  // #define PETSC_CUPMBLAS_FP_INPUT_TYPE  D
  // #define PETSC_CUPMBLAS_FP_RETURN_TYPE z
  // PETSC_CUPMBLAS_BUILD_BLAS_FUNCTION_ALIAS_MODIFIED(nrm2) -> Dznrm2
  #define PETSC_CUPMBLAS_BUILD_BLAS_FUNCTION_ALIAS_MODIFIED(func) PetscConcat(PetscConcat(PETSC_CUPMBLAS_FP_INPUT_TYPE, PETSC_CUPMBLAS_FP_RETURN_TYPE), func)

  // PETSC_CUPMBLAS_BUILD_BLAS_FUNCTION_ALIAS_IFPTYPE() - Helper macro to build Iamax and Iamin
  // because they are both extra special
  //
  // input param:
  // func - base suffix of the blas function, either amax or amin
  //
  // notes:
  // The macro name literally stands for "I" ## "floating point type" because shockingly enough,
  // that's what it does.
  //
  // requires PETSC_CUPMBLAS_FP_TYPE_L to be defined as the lower-case blas floating point input type
  // letter ("s" for complex single, "z" for complex double, "s" for real single, and "d" for
  // real double).
  //
  // example usage:
  // #define PETSC_CUPMBLAS_FP_TYPE_L s
  // PETSC_CUPMBLAS_BUILD_BLAS_FUNCTION_ALIAS_IFPTYPE(amax) -> Isamax
  //
  // #define PETSC_CUPMBLAS_FP_TYPE_L z
  // PETSC_CUPMBLAS_BUILD_BLAS_FUNCTION_ALIAS_IFPTYPE(amin) -> Izamin
  #define PETSC_CUPMBLAS_BUILD_BLAS_FUNCTION_ALIAS_IFPTYPE(func) PetscConcat(I, PetscConcat(PETSC_CUPMBLAS_FP_TYPE_L, func))

  // PETSC_CUPMBLAS_BUILD_BLAS_FUNCTION_ALIAS_STANDARD() - Helper macro to build a "standard"
  // blas function name
  //
  // input param:
  // func - base suffix of the blas function, e.g. axpy, scal
  //
  // notes:
  // requires PETSC_CUPMBLAS_FP_TYPE to be defined as the blas floating-point letter ("C" for
  // complex single, "Z" for complex double, "S" for real single, "D" for real double).
  //
  // example usage:
  // #define PETSC_CUPMBLAS_FP_TYPE S
  // PETSC_CUPMBLAS_BUILD_BLAS_FUNCTION_ALIAS_STANDARD(axpy) -> Saxpy
  //
  // #define PETSC_CUPMBLAS_FP_TYPE Z
  // PETSC_CUPMBLAS_BUILD_BLAS_FUNCTION_ALIAS_STANDARD(axpy) -> Zaxpy
  #define PETSC_CUPMBLAS_BUILD_BLAS_FUNCTION_ALIAS_STANDARD(func) PetscConcat(PETSC_CUPMBLAS_FP_TYPE, func)

  // PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION_EXACT() - In case CUDA/HIP don't agree with our suffix
  // one can provide both here
  //
  // input params:
  // MACRO_SUFFIX - suffix to one of the above blas function builder macros, e.g. STANDARD or
  // IFPTYPE
  // our_suffix   - the suffix of the alias function
  // their_suffix - the suffix of the funciton being aliased
  //
  // notes:
  // requires PETSC_CUPMBLAS_PREFIX to be defined as the specific CUDA/HIP blas function
  // prefix. requires any other specific definitions required by the specific builder macro to
  // also be defined. See PETSC_CUPM_ALIAS_FUNCTION_EXACT() for the exact expansion of the
  // function alias.
  //
  // example usage:
  // #define PETSC_CUPMBLAS_PREFIX  cublas
  // #define PETSC_CUPMBLAS_FP_TYPE C
  // PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION_EXACT(STANDARD,dot,dotc) ->
  // template <typename... T>
  // static constexpr auto cupmBlasXdot(T&&... args) *noexcept and returntype detection*
  // {
  //   return cublasCdotc(std::forward<T>(args)...);
  // }
  #define PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION_EXACT(MACRO_SUFFIX, our_suffix, their_suffix) \
    PETSC_CUPM_ALIAS_FUNCTION(PetscConcat(cupmBlasX, our_suffix), PetscConcat(PETSC_CUPMBLAS_PREFIX, PetscConcat(PETSC_CUPMBLAS_BUILD_BLAS_FUNCTION_ALIAS_, MACRO_SUFFIX)(their_suffix)))

  // PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION() - Alias a CUDA/HIP blas function
  //
  // input params:
  // MACRO_SUFFIX - suffix to one of the above blas function builder macros, e.g. STANDARD or
  // IFPTYPE
  // suffix       - the common suffix between CUDA and HIP of the alias function
  //
  // notes:
  // see PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(), this macro just calls that one with "suffix" as
  // "our_prefix" and "their_prefix"
  #define PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(MACRO_SUFFIX, suffix) PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION_EXACT(MACRO_SUFFIX, suffix, suffix)

  // PETSC_CUPMBLAS_ALIAS_FUNCTION() - Alias a CUDA/HIP library function
  //
  // input params:
  // suffix - the common suffix between CUDA and HIP of the alias function
  //
  // notes:
  // requires PETSC_CUPMBLAS_PREFIX to be defined as the specific CUDA/HIP blas library
  // prefix. see PETSC_CUPMM_ALIAS_FUNCTION_EXACT() for the precise expansion of this macro.
  //
  // example usage:
  // #define PETSC_CUPMBLAS_PREFIX hipblas
  // PETSC_CUPMBLAS_ALIAS_FUNCTION(Create) ->
  // template <typename... T>
  // static constexpr auto cupmBlasCreate(T&&... args) *noexcept and returntype detection*
  // {
  //   return hipblasCreate(std::forward<T>(args)...);
  // }
  #define PETSC_CUPMBLAS_ALIAS_FUNCTION(suffix) PETSC_CUPM_ALIAS_FUNCTION(PetscConcat(cupmBlas, suffix), PetscConcat(PETSC_CUPMBLAS_PREFIX, suffix))

template <DeviceType T>
struct BlasInterfaceBase : Interface<T> {
  PETSC_NODISCARD static constexpr const char *cupmBlasName() noexcept { return T == DeviceType::CUDA ? "cuBLAS" : "hipBLAS"; }
};

  #define PETSC_CUPMBLAS_BASE_CLASS_HEADER(DEV_TYPE) \
    using base_type = ::Petsc::device::cupm::impl::BlasInterfaceBase<DEV_TYPE>; \
    using base_type::cupmBlasName; \
    PETSC_CUPM_ALIAS_FUNCTION(cupmBlasGetErrorName, PetscConcat(PetscConcat(Petsc, PETSC_CUPMBLAS_PREFIX_U), GetErrorName)) \
    PETSC_CUPM_INHERIT_INTERFACE_TYPEDEFS_USING(interface_type, DEV_TYPE)

template <DeviceType>
struct BlasInterfaceImpl;

  #if PetscDefined(HAVE_CUDA)
    #define PETSC_CUPMBLAS_PREFIX         cublas
    #define PETSC_CUPMBLAS_PREFIX_U       CUBLAS
    #define PETSC_CUPMBLAS_FP_TYPE        PETSC_CUPMBLAS_FP_TYPE_U
    #define PETSC_CUPMBLAS_FP_INPUT_TYPE  PETSC_CUPMBLAS_FP_INPUT_TYPE_U
    #define PETSC_CUPMBLAS_FP_RETURN_TYPE PETSC_CUPMBLAS_FP_RETURN_TYPE_L
template <>
struct BlasInterfaceImpl<DeviceType::CUDA> : BlasInterfaceBase<DeviceType::CUDA> {
  PETSC_CUPMBLAS_BASE_CLASS_HEADER(DeviceType::CUDA);

  // typedefs
  using cupmBlasHandle_t      = cublasHandle_t;
  using cupmBlasError_t       = cublasStatus_t;
  using cupmBlasInt_t         = int;
  using cupmSolverHandle_t    = cusolverDnHandle_t;
  using cupmSolverError_t     = cusolverStatus_t;
  using cupmBlasPointerMode_t = cublasPointerMode_t;

  // values
  static const auto CUPMBLAS_STATUS_SUCCESS         = CUBLAS_STATUS_SUCCESS;
  static const auto CUPMBLAS_STATUS_NOT_INITIALIZED = CUBLAS_STATUS_NOT_INITIALIZED;
  static const auto CUPMBLAS_STATUS_ALLOC_FAILED    = CUBLAS_STATUS_ALLOC_FAILED;
  static const auto CUPMBLAS_POINTER_MODE_HOST      = CUBLAS_POINTER_MODE_HOST;
  static const auto CUPMBLAS_POINTER_MODE_DEVICE    = CUBLAS_POINTER_MODE_DEVICE;

  // utility functions
  PETSC_CUPMBLAS_ALIAS_FUNCTION(Create)
  PETSC_CUPMBLAS_ALIAS_FUNCTION(Destroy)
  PETSC_CUPMBLAS_ALIAS_FUNCTION(GetStream)
  PETSC_CUPMBLAS_ALIAS_FUNCTION(SetStream)
  PETSC_CUPMBLAS_ALIAS_FUNCTION(GetPointerMode)
  PETSC_CUPMBLAS_ALIAS_FUNCTION(SetPointerMode)

  // level 1 BLAS
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(STANDARD, axpy)
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(STANDARD, scal)
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION_EXACT(STANDARD, dot, PetscIfPetscDefined(USE_COMPLEX, dotc, dot))
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION_EXACT(STANDARD, dotu, PetscIfPetscDefined(USE_COMPLEX, dotu, dot))
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(STANDARD, swap)
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(MODIFIED, nrm2)
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(IFPTYPE, amax)
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(MODIFIED, asum)

  // level 2 BLAS
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(STANDARD, gemv)

  // level 3 BLAS
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(STANDARD, gemm)

  // BLAS extensions
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(STANDARD, geam)

  PETSC_NODISCARD static PetscErrorCode InitializeHandle(cupmSolverHandle_t &handle) noexcept
  {
    PetscFunctionBegin;
    if (handle) PetscFunctionReturn(0);
    for (auto i = 0; i < 3; ++i) {
      const auto cerr = cusolverDnCreate(&handle);
      if (PetscLikely(cerr == CUSOLVER_STATUS_SUCCESS)) break;
      if ((cerr != CUSOLVER_STATUS_NOT_INITIALIZED) && (cerr != CUSOLVER_STATUS_ALLOC_FAILED)) PetscCallCUSOLVER(cerr);
      if (i < 2) {
        PetscCall(PetscSleep(3));
        continue;
      }
      PetscCheck(cerr == CUSOLVER_STATUS_SUCCESS, PETSC_COMM_SELF, PETSC_ERR_GPU_RESOURCE, "Unable to initialize cuSolverDn");
    }
    PetscFunctionReturn(0);
  }

  PETSC_NODISCARD static PetscErrorCode SetHandleStream(const cupmSolverHandle_t &handle, const cupmStream_t &stream) noexcept
  {
    cupmStream_t cupmStream;

    PetscFunctionBegin;
    PetscCallCUSOLVER(cusolverDnGetStream(handle, &cupmStream));
    if (cupmStream != stream) PetscCallCUSOLVER(cusolverDnSetStream(handle, stream));
    PetscFunctionReturn(0);
  }

  PETSC_NODISCARD static PetscErrorCode DestroyHandle(cupmSolverHandle_t &handle) noexcept
  {
    PetscFunctionBegin;
    if (handle) {
      PetscCallCUSOLVER(cusolverDnDestroy(handle));
      handle = nullptr;
    }
    PetscFunctionReturn(0);
  }
};
    #undef PETSC_CUPMBLAS_PREFIX
    #undef PETSC_CUPMBLAS_PREFIX_U
    #undef PETSC_CUPMBLAS_FP_TYPE
    #undef PETSC_CUPMBLAS_FP_INPUT_TYPE
    #undef PETSC_CUPMBLAS_FP_RETURN_TYPE
  #endif // PetscDefined(HAVE_CUDA)

  #if PetscDefined(HAVE_HIP)
    #define PETSC_CUPMBLAS_PREFIX         hipblas
    #define PETSC_CUPMBLAS_PREFIX_U       HIPBLAS
    #define PETSC_CUPMBLAS_FP_TYPE        PETSC_CUPMBLAS_FP_TYPE_U
    #define PETSC_CUPMBLAS_FP_INPUT_TYPE  PETSC_CUPMBLAS_FP_INPUT_TYPE_U
    #define PETSC_CUPMBLAS_FP_RETURN_TYPE PETSC_CUPMBLAS_FP_RETURN_TYPE_L
template <>
struct BlasInterfaceImpl<DeviceType::HIP> : BlasInterfaceBase<DeviceType::HIP> {
  PETSC_CUPMBLAS_BASE_CLASS_HEADER(DeviceType::HIP);

  // typedefs
  using cupmBlasHandle_t      = hipblasHandle_t;
  using cupmBlasError_t       = hipblasStatus_t;
  using cupmBlasInt_t         = int; // rocblas will have its own
  using cupmSolverHandle_t    = hipsolverHandle_t;
  using cupmSolverError_t     = hipsolverStatus_t;
  using cupmBlasPointerMode_t = hipblasPointerMode_t;

  // values
  static const auto CUPMBLAS_STATUS_SUCCESS         = HIPBLAS_STATUS_SUCCESS;
  static const auto CUPMBLAS_STATUS_NOT_INITIALIZED = HIPBLAS_STATUS_NOT_INITIALIZED;
  static const auto CUPMBLAS_STATUS_ALLOC_FAILED    = HIPBLAS_STATUS_ALLOC_FAILED;
  static const auto CUPMBLAS_POINTER_MODE_HOST      = HIPBLAS_POINTER_MODE_HOST;
  static const auto CUPMBLAS_POINTER_MODE_DEVICE    = HIPBLAS_POINTER_MODE_DEVICE;

  // utility functions
  PETSC_CUPMBLAS_ALIAS_FUNCTION(Create)
  PETSC_CUPMBLAS_ALIAS_FUNCTION(Destroy)
  PETSC_CUPMBLAS_ALIAS_FUNCTION(GetStream)
  PETSC_CUPMBLAS_ALIAS_FUNCTION(SetStream)
  PETSC_CUPMBLAS_ALIAS_FUNCTION(GetPointerMode)
  PETSC_CUPMBLAS_ALIAS_FUNCTION(SetPointerMode)

  // level 1 BLAS
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(STANDARD, axpy)
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(STANDARD, scal)
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION_EXACT(STANDARD, dot, PetscIfPetscDefined(USE_COMPLEX, dotc, dot))
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION_EXACT(STANDARD, dotu, PetscIfPetscDefined(USE_COMPLEX, dotu, dot))
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(STANDARD, swap)
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(MODIFIED, nrm2)
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(IFPTYPE, amax)
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(MODIFIED, asum)

  // level 2 BLAS
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(STANDARD, gemv)

  // level 3 BLAS
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(STANDARD, gemm)

  // BLAS extensions
  PETSC_CUPMBLAS_ALIAS_BLAS_FUNCTION(STANDARD, geam)

  PETSC_NODISCARD static PetscErrorCode InitializeHandle(cupmSolverHandle_t &handle) noexcept
  {
    PetscFunctionBegin;
    if (!handle) PetscCallHIPSOLVER(hipsolverCreate(&handle));
    PetscFunctionReturn(0);
  }

  PETSC_NODISCARD static PetscErrorCode SetHandleStream(cupmSolverHandle_t handle, cupmStream_t stream) noexcept
  {
    PetscFunctionBegin;
    PetscCallHIPSOLVER(hipsolverSetStream(handle, stream));
    PetscFunctionReturn(0);
  }

  PETSC_NODISCARD static PetscErrorCode DestroyHandle(cupmSolverHandle_t &handle) noexcept
  {
    PetscFunctionBegin;
    if (handle) {
      PetscCallHIPSOLVER(hipsolverDestroy(handle));
      handle = nullptr;
    }
    PetscFunctionReturn(0);
  }
};
    #undef PETSC_CUPMBLAS_PREFIX
    #undef PETSC_CUPMBLAS_PREFIX_U
    #undef PETSC_CUPMBLAS_FP_TYPE
    #undef PETSC_CUPMBLAS_FP_INPUT_TYPE
    #undef PETSC_CUPMBLAS_FP_RETURN_TYPE
  #endif // PetscDefined(HAVE_HIP)

  #undef PETSC_CUPMBLAS_BASE_CLASS_HEADER

  #define PETSC_CUPMBLAS_IMPL_CLASS_HEADER(base_name, T) \
    PETSC_CUPM_INHERIT_INTERFACE_TYPEDEFS_USING(cupmInterface_t, T); \
    using base_name = ::Petsc::device::cupm::impl::BlasInterfaceImpl<T>; \
    /* introspection */ \
    using base_name::cupmBlasName; \
    using base_name::cupmBlasGetErrorName; \
    /* types */ \
    using cupmBlasHandle_t      = typename base_name::cupmBlasHandle_t; \
    using cupmBlasError_t       = typename base_name::cupmBlasError_t; \
    using cupmBlasInt_t         = typename base_name::cupmBlasInt_t; \
    using cupmSolverHandle_t    = typename base_name::cupmSolverHandle_t; \
    using cupmSolverError_t     = typename base_name::cupmSolverError_t; \
    using cupmBlasPointerMode_t = typename base_name::cupmBlasPointerMode_t; \
    /* values */ \
    using base_name::CUPMBLAS_STATUS_SUCCESS; \
    using base_name::CUPMBLAS_STATUS_NOT_INITIALIZED; \
    using base_name::CUPMBLAS_STATUS_ALLOC_FAILED; \
    using base_name::CUPMBLAS_POINTER_MODE_HOST; \
    using base_name::CUPMBLAS_POINTER_MODE_DEVICE; \
    /* utility functions */ \
    using base_name::cupmBlasCreate; \
    using base_name::cupmBlasDestroy; \
    using base_name::cupmBlasGetStream; \
    using base_name::cupmBlasSetStream; \
    using base_name::cupmBlasGetPointerMode; \
    using base_name::cupmBlasSetPointerMode; \
    /* level 1 BLAS */ \
    using base_name::cupmBlasXaxpy; \
    using base_name::cupmBlasXscal; \
    using base_name::cupmBlasXdot; \
    using base_name::cupmBlasXdotu; \
    using base_name::cupmBlasXswap; \
    using base_name::cupmBlasXnrm2; \
    using base_name::cupmBlasXamax; \
    using base_name::cupmBlasXasum; \
    /* level 2 BLAS */ \
    using base_name::cupmBlasXgemv; \
    /* level 3 BLAS */ \
    using base_name::cupmBlasXgemm; \
    /* BLAS extensions */ \
    using base_name::cupmBlasXgeam

// The actual interface class
template <DeviceType T>
struct BlasInterface : BlasInterfaceImpl<T> {
  PETSC_CUPMBLAS_IMPL_CLASS_HEADER(blasinterface_type, T);

  PETSC_NODISCARD static PetscErrorCode PetscCUPMBlasSetPointerModeFromPointer(cupmBlasHandle_t handle, const void *ptr) noexcept
  {
    auto mtype = PETSC_MEMTYPE_HOST;

    PetscFunctionBegin;
    PetscCall(PetscCUPMGetMemType(ptr, &mtype));
    PetscCallCUPMBLAS(cupmBlasSetPointerMode(handle, PetscMemTypeDevice(mtype) ? CUPMBLAS_POINTER_MODE_DEVICE : CUPMBLAS_POINTER_MODE_HOST));
    PetscFunctionReturn(0);
  }
};

  #define PETSC_CUPMBLAS_INHERIT_INTERFACE_TYPEDEFS_USING(base_name, T) \
    PETSC_CUPMBLAS_IMPL_CLASS_HEADER(PetscConcat(base_name, _impl), T); \
    using base_name = ::Petsc::device::cupm::impl::BlasInterface<T>; \
    using base_name::PetscCUPMBlasSetPointerModeFromPointer

  #if PetscDefined(HAVE_CUDA)
extern template struct BlasInterface<DeviceType::CUDA>;
  #endif

  #if PetscDefined(HAVE_HIP)
extern template struct BlasInterface<DeviceType::HIP>;
  #endif

} // namespace impl

} // namespace cupm

} // namespace device

} // namespace Petsc

#endif // defined(__cplusplus)

#endif // PETSCCUPMBLASINTERFACE_HPP
