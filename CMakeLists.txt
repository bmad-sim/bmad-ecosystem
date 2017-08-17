cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set (LIBNAME tao)

set (CREATE_SHARED true)

SET (PLOT_LINK_LIBS $ENV{PLOT_LINK_LIBS})

## Profiler flags
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# set (FFLAGS "-pg")  # For profiling
# set (CFLAGS "-pg")

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

IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  SET (SHARED_DEPS
    xsif
    xrlf03
    xrl
    sim_utils
    recipes_f-90_LEPP 
    ${PLOT_LINK_LIBS}
    forest 
    fgsl
    gsl
    gslcblas
    ${SHARED_DEPS}
    X11
    readline
  )
ENDIF ()

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
