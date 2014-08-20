cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(LIBNAME sim_utils)

set (CREATE_SHARED true)

set(INC_DIRS 
)

set(SRC_DIRS
  file 
  interfaces 
  io 
  math 
  matrix 
  misc 
  plot 
  string
  geodesic_lm
)

set(TEST_EXE_SPECS
  cmake.test
)

# set (FFLAGS "-pg")  # For profiling
# set (CFLAGS "-pg")

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





