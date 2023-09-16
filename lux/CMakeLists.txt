cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

set (LIBNAME lux)

set (SRC_DIRS
  code
)

set(EXE_SPECS cmake.lux)

if ($ENV{ACC_ENABLE_MPI})
  message ('Building MPI version...')
  list (APPEND EXE_SPECS cmake.lux_mpi)
endif ()

# For profiling...
# set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# set (FFLAGS "-pg") 
# set (CFLAGS "-pg")

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
