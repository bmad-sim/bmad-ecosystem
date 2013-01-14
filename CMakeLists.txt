cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set (LIBNAME tao)

set(INC_DIRS 
)

set(SRC_DIRS
  code 
  hook
)

set (EXE_SPECS cmake.tao)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





