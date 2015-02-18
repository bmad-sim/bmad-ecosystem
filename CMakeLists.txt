cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

if ($ENV{ACC_ENABLE_MPI})
  message ('Building MPI version...')
  set(EXE_SPECS cmake.lux_mpi)
else ()
  set(EXE_SPECS cmake.lux)
endif ()

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
