
set(LIBNAME bmad)

cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set (CREATE_SHARED true)

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
  low_level
)

# For profiling...
# set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# set (FFLAGS "-pg") 
# set (CFLAGS "-pg")

# For SHARED support
SET (FFLAGS "-fPIC")
SET (LINK_FLAGS "-fPIC")

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





