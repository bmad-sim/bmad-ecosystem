cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

SET (LIBNAME sim_utils)

SET (CREATE_SHARED true)

SET (PLOT_LINK_LIBS $ENV{PLOT_LINK_LIBS})

# If this is a conda build, fgsl is contained
# one level deeper in the include folder
IF ($ENV{ACC_CONDA_BUILD})
  SET (INC_DIRS
    $ENV{ACC_CONDA_PATH}/include/
    $ENV{ACC_CONDA_PATH}/include/fgsl/
  )
ELSE ()
  SET (INC_DIRS
  )
ENDIF ()

SET (SRC_DIRS
  file 
  interfaces 
  io 
  math 
  matrix 
  misc
  optimizers
  plot 
  string
  geodesic_lm
  special_functions
  )

SET (SHARED_DEPS
  ${PLOT_LINK_LIBS}
  $ENV{ACC_BMAD_LINK_LIBS}
  X11
  readline
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
