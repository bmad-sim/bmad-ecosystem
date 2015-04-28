
set(LIBNAME bmad)

cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(INC_DIRS
)

set(SRC_DIRS 
  code 
  custom
  interface
  modules
  multiparticle
  photon
  parsing
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





