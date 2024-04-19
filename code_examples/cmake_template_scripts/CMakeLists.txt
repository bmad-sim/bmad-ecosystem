# This is a minimal example of a CMakeLists.txt file to make an executable.
# This script works in conjunction with the script cmake.test that 
# is included in the same directory as this file.

cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

set(EXE_SPECS cmake.test)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
