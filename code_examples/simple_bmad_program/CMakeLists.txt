cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

set(EXE_SPECS 
  cmake.simple_bmad_program
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
