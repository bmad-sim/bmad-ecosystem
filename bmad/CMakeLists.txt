cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

SET (LIBNAME bmad)

SET (CREATE_SHARED true)

SET (PLOT_LINK_LIBS $ENV{PLOT_LINK_LIBS})

# If this is a conda build, xraylib is contained
# one level deeper in the include folder

IF ($ENV{ACC_CONDA_BUILD})
  SET (INC_DIRS
    $ENV{ACC_CONDA_PATH}/include
    $ENV{ACC_CONDA_PATH}/include/xraylib/
    $ENV{ACC_CONDA_PATH}/include/fgsl/
  )
ELSE ()
  SET (INC_DIRS
  )
ENDIF ()

# ---- GPU-accelerated tracking for drift and quadrupole ----
# Enable by setting ACC_ENABLE_GPU_TRACKING=Y in the environment.
# Requires CUDA Toolkit (nvcc for kernel compilation, libcudart for runtime).
IF ("$ENV{ACC_ENABLE_GPU_TRACKING}" STREQUAL "Y")
  enable_language(CUDA)
  find_package(CUDAToolkit QUIET)
  IF (CUDAToolkit_FOUND AND CMAKE_CUDA_COMPILER)
    MESSAGE(STATUS "GPU tracking enabled — CUDA ${CUDAToolkit_VERSION}, nvcc: ${CMAKE_CUDA_COMPILER}")
    LIST(APPEND INC_DIRS ${CUDAToolkit_INCLUDE_DIRS})
    add_definitions(-DUSE_GPU_TRACKING)
    # Set CUDA architecture — override via ACC_GPU_CUDA_ARCH env var (default: native)
    IF (DEFINED ENV{ACC_GPU_CUDA_ARCH})
      SET(CMAKE_CUDA_ARCHITECTURES $ENV{ACC_GPU_CUDA_ARCH})
    ELSE ()
      SET(CMAKE_CUDA_ARCHITECTURES native)
    ENDIF ()
    MESSAGE(STATUS "GPU tracking CUDA architectures: ${CMAKE_CUDA_ARCHITECTURES}")
    # Build CUDA kernels as a separate static library
    add_library(gpu_tracking_kernels STATIC
      ${CMAKE_CURRENT_SOURCE_DIR}/low_level/gpu_tracking_kernels.cu)
    target_compile_definitions(gpu_tracking_kernels PRIVATE USE_GPU_TRACKING)
    set_target_properties(gpu_tracking_kernels PROPERTIES
      POSITION_INDEPENDENT_CODE ON)
    SET (GPU_TRACKING_LINK_LIBS gpu_tracking_kernels CUDA::cudart)
  ELSE ()
    MESSAGE(WARNING "ACC_ENABLE_GPU_TRACKING=Y but CUDA/nvcc not found — GPU tracking disabled")
  ENDIF ()
ENDIF ()

SET (SRC_DIRS 
  code 
  ### custom  # Using function pointers now
  geometry
  interface
  modules
  multiparticle
  output
  ptc
  photon
  parsing
  low_level
  hdf5
  space_charge
  spin
  )

SET (SHARED_DEPS
  sim_utils
  $ENV{ACC_BMAD_LINK_LIBS}
  ${GPU_TRACKING_LINK_LIBS}
  X11
  readline
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
