cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

if ($ENV{ACC_ENABLE_MPI})
  message ('Building MPI version...')
  set(EXE_SPECS cmake.lux_mpi)
else ()
  set(EXE_SPECS cmake.lux)
endif ()

# For profiling...
# set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# set (FFLAGS "-pg") 
# set (CFLAGS "-pg")

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
