cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set (EXE_SPECS 
  cmake.abs_time_test
  cmake.autoscale_test
  cmake.coord_test
  cmake.bookkeeper_test
  cmake.geometry_test
  cmake.girder_test
  cmake.mat6_calc_method_test
  cmake.multipass_test
  cmake.parse_test
  cmake.patch_test
  cmake.photon_test
  cmake.ptc_test
  cmake.reverse_test
  cmake.sim_utils_test
  cmake.slice_test
  cmake.superimpose_test
  cmake.synrad3d_test
  cmake.time_runge_kutta_test
  cmake.tracking_method_test
  cmake.twiss_track_test
  cmake.wake_test
  cmake.wall3d_test
  cmake.xraylib_test
)

set(SRC_DIRS
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





