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

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
