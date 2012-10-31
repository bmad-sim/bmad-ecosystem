
set(LIBNAME bmad)

cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(INC_DIRS
  ../include
  $ENV{ACC_PKG}/fgsl
)

set(SRC_DIRS 
  code 
  custom
  modules
  synrad
  photon
  parsing
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





