cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(INC_DIRS 
  include
)

set(SRC_DIRS
  interface_test
  util
)

set(EXE_SPECS interface_test.cmake)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





