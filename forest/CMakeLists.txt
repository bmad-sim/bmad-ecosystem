cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

SET (LIBNAME forest)

SET (CREATE_SHARED true)

SET (INC_DIRS 
  include
  )

SET (SRC_DIRS
  code
  )

## Profiler flags - Note: This  is not compatible with the -fPIC flag, so -pg can't be set for SHARED build.
# set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# set (FFLAGS "-pg") 
# set (CFLAGS "-pg")

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





