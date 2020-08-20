cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

set (EXE_SPECS 
  cmake.abs_time_test
  cmake.aperture_test
  cmake.autoscale_test
  cmake.bbu_test
  cmake.closed_orbit_test
  cmake.coord_test
  cmake.csr_and_space_charge
  cmake.bookkeeper_test
  cmake.em_field_test
  cmake.field_test
  cmake.geometry_test
  cmake.girder_test
  cmake.hdf5_test
  cmake.mat6_calc_method_test
  cmake.mode_flip_test
  cmake.mode3_test
  cmake.multipass_test
  cmake.multipole_test
  cmake.object_test
  cmake.parse_test
  cmake.particle_species_test
  cmake.patch_test
  cmake.photon_test
  cmake.ptc_test
  cmake.reverse_test
  cmake.sad_test
  cmake.sim_utils_test
  cmake.slice_test
  cmake.spin_test
  cmake.superimpose_test
  cmake.synrad3d_test
  cmake.taylor_test
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





