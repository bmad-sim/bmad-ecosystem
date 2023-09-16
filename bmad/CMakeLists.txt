cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

SET (LIBNAME bmad)

SET (CREATE_SHARED true)

SET (PLOT_LINK_LIBS $ENV{PLOT_LINK_LIBS})

# If this is a conda build, xraylib is contained
# one level deeper in the include folder
IF ($ENV{ACC_CONDA_BUILD})
  SET (INC_DIRS
    $ENV{ACC_CONDA_PATH}/include/xraylib/
    $ENV{ACC_CONDA_PATH}/include/fgsl/
  )
ELSE ()
  SET (INC_DIRS
  )
ENDIF ()

SET (SRC_DIRS 
  code 
  custom
  geometry
  interface
  modules
  multiparticle
  ptc
  photon
  parsing
  low_level
  hdf5
  space_charge
  spin
  )

# For profiling... - Note: This is not compatible with the -fPIC flag, so -pg can't be set for SHARED build.
# SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# SET (FFLAGS "-pg") 
# SET (CFLAGS "-pg")

# For SHARED support - Note: should now be set automatically by Master.cmake
# SET (FFLAGS "-fPIC")
# SET (LINK_FLAGS "-fPIC")

SET (SHARED_DEPS
  sim_utils
  $ENV{ACC_BMAD_LINK_LIBS}
  X11
  readline
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
