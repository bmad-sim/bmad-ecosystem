
set (LIBNAME tao)
cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(INC_DIRS 
  ../include
  include
)

set(SRC_DIRS
  code 
  hook
)

set (EXE_SPECS tao.cmake)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





