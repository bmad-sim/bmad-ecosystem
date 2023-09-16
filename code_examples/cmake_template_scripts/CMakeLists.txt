# This is a minimal example of a CMakeLists.txt file to make an executable.
# This script works in conjunction with the script cmake.test that 
# is included in the same directory as this file.

cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

set(EXE_SPECS cmake.test)

## Profiler flags - Note: This	is not compatible with the -fPIC flag, so -pg can't be set for SHARED build.
# SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# SET (FFLAGS "-pg")  # For profiling
# SET (CFLAGS "-pg")

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
