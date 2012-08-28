cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set (EXE_SPECS 
  ptc_test.cmake
  time_runge_kutta_test.cmake
  synrad3d_test.cmake
  twiss_track_test.cmake
  slice_test.cmake
  abs_time_test.cmake
)

set(SRC_DIRS
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





