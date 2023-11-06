cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

SET (LIBNAME tao)

SET (CREATE_SHARED true)

SET (PLOT_LINK_LIBS $ENV{PLOT_LINK_LIBS})

## Profiler flags - Note: This	is not compatible with the -fPIC flag, so -pg can't be set for SHARED build.
# SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# SET (FFLAGS "-pg")  # For profiling
# SET (CFLAGS "-pg")

SET (SRC_DIRS
  code 
  ## hook # Now using function pointers
  version
  )

SET (EXE_SPECS
  cmake.tao
  )

SET (SHARED_DEPS
  bmad
  sim_utils
  ${PLOT_LINK_LIBS}
  $ENV{ACC_BMAD_LINK_LIBS}
  gomp
  X11
  readline
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
