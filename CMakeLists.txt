
set(LIBNAME bmad)

cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set (CREATE_SHARED true)

set(INC_DIRS
)

set(SRC_DIRS 
  code 
  custom
  interface
  modules
  multiparticle
  photon
  parsing
  low_level
)

# For profiling...
# set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# set (FFLAGS "-pg") 
# set (CFLAGS "-pg")

# For SHARED support
SET (FFLAGS "-fPIC")
SET (LINK_FLAGS "-fPIC")

IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  SET (SHARED_DEPS
    bmad 
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
    lapack95
    lapack
    blas
    X11
    readline
  )
ENDIF ()

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





