cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set (EXE_SPECS 
  cmake.girder_test
  cmake.parse_test
  cmake.reverse_test
  cmake.multipass_test
  cmake.sim_utils_test
  cmake.ptc_test
  cmake.time_runge_kutta_test
  cmake.synrad3d_test
  cmake.twiss_track_test
  cmake.slice_test
  cmake.abs_time_test
  cmake.tracking_method_test
  cmake.mat6_calc_method_test
)

set(SRC_DIRS
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





