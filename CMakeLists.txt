cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set (LIBNAME tao)

## set (CREATE_SHARED true)
## Profiler flags
## set (FFLAGS  "-pg")

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

set (FFLAGS "-pg")  # For profiling
set (CFLAGS "-pg")

set(SRC_DIRS
  code 
  hook
)

set (EXE_SPECS
  cmake.tao
)

set(SHARED_DEPS
  bmad
  lapack95
  lapack
  blas
  )

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





