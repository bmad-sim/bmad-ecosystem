cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

set (EXE_SPECS 
  cmake_files/cmake.abs_time_test
  cmake_files/cmake.analysis_test
  cmake_files/cmake.aperture_test
  cmake_files/cmake.autoscale_test
  cmake_files/cmake.backwards_time_track_test
  cmake_files/cmake.beam_test
  cmake_files/cmake.bbu_test
  cmake_files/cmake.cesr_test
  cmake_files/cmake.closed_orbit_test
  cmake_files/cmake.coord_test
  cmake_files/cmake.csr_and_space_charge
  cmake_files/cmake.bookkeeper_test
  cmake_files/cmake.cathode_sc_test
  cmake_files/cmake.em_field_test
  cmake_files/cmake.field_test
  cmake_files/cmake.fringe_test
  cmake_files/cmake.geometry_test
  cmake_files/cmake.girder_test
  cmake_files/cmake.hdf5_test
  cmake_files/cmake.mad_to_bmad_test
  cmake_files/cmake.mat6_calc_method_test
  cmake_files/cmake.match_test
  cmake_files/cmake.mode_flip_test
  cmake_files/cmake.mode3_test
  cmake_files/cmake.multipass_test
  cmake_files/cmake.multipole_test
  cmake_files/cmake.nonlin_test
  cmake_files/cmake.object_test
  cmake_files/cmake.parse_test
  cmake_files/cmake.particle_species_test
  cmake_files/cmake.patch_test
  cmake_files/cmake.photon_test
  cmake_files/cmake.ptc_test
  cmake_files/cmake.radiation_test
  cmake_files/cmake.reverse_test
  cmake_files/cmake.sad_test
  cmake_files/cmake.sim_utils_test
  cmake_files/cmake.slice_test
  cmake_files/cmake.space_charge_test
  cmake_files/cmake.spin_map_test
  cmake_files/cmake.spin_general_test
  cmake_files/cmake.superimpose_test
  cmake_files/cmake.synrad3d_test
  cmake_files/cmake.taylor_test
  cmake_files/cmake.time_runge_kutta_test
  cmake_files/cmake.tracking_method_test
  cmake_files/cmake.twiss_test
  cmake_files/cmake.wake_test
  cmake_files/cmake.wall3d_test
  cmake_files/cmake.write_bmad_test
  cmake_files/cmake.xraylib_test
)

set(SRC_DIRS
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





