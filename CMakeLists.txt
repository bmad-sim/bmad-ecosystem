cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(LIBNAME sim_utils)

set(INC_DIRS 
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





