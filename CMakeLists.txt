cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set (LIBNAME tao)

set (CREATE_SHARED true)

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
  IF ($ENV{ACC_PLOT_PACKAGE} MATCHES "pgplot")
    SET (PLOT_LINK_LIBS "pgplot")
  ELSEIF ($ENV{ACC_PLOT_PACKAGE} MATCHES "plplot")
    SET (PLOT_LINK_LIBS plplot plplotcxx plplotf95 csirocsa qsastime)
  ELSEIF ($ENV{ACC_PLOT_PACKAGE} MATCHES "none")
    SET (PLOT_LINK_LIBS)
  ENDIF ()

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
