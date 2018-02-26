set(LIBNAME forest)
cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(INC_DIRS 
  include
)

set(SRC_DIRS
  code
)

# For profiling...
# set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# set (FFLAGS "-pg") 
# set (CFLAGS "-pg")

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)






