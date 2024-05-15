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
  X11
  readline
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
