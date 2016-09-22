cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set (LIBNAME cpp_bmad_interface)

set(INC_DIRS 
  include
)

set(SRC_DIRS
  code
)

set(EXE_SPECS 
  cmake.cpp_bmad_interface_test
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)




