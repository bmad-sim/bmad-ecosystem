cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(INC_DIRS 
  ../include
  include
)

set(SRC_DIRS
)

set (EXE_SPECS ptc_test.cmake)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





