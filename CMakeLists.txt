
set(LIBNAME bmad)

cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(INC_DIRS
  ../include 
  include
)

set(SRC_DIRS 
  code 
  custom
  modules
  cpp_interface
  synrad
  photon
  parsing
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





