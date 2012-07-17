
set(LIBNAME sim_utils)

cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(INC_DIRS 
  ../include
)

set(SRC_DIRS
  file 
  interfaces 
  io 
  math 
  matrix 
  misc 
  old_utils 
  plot 
  string
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





