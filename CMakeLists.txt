cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

SET (LIBNAME sim_utils)

SET (CREATE_SHARED true)

SET (PLOT_LINK_LIBS $ENV{PLOT_LINK_LIBS})

SET (INC_DIRS 
  )

SET (SRC_DIRS
  file 
  code
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

## Profiler flags - Note: This  is not compatible with the -fPIC flag, so -pg can't be set for the SHARED build.
# SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# SET (FFLAGS "-pg")  
# SET (CFLAGS "-pg")

SET (SHARED_DEPS
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

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
