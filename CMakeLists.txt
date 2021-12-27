set (LIBNAME bsim)
cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

set(INC_DIRS 
)

set (SRC_DIRS
  code
  code_synrad
  code_synrad3d
  synrad3d/custom
  modules
)

set (EXE_SPECS
  #cmake.ansga2
  #cmake.aspea2
  cmake.dark_current_tracker
  cmake.envelope_ibs
  cmake.srdt_lsq_soln
  cmake.frequency_map
  cmake.bbu
  cmake.dynamic_aperture
  cmake.ibs_linac
  cmake.ibs_ring
  cmake.spin_stroboscope
  cmake.synrad
  cmake.synrad3d
  cmake.synrad_aperture_to_wall
  cmake.tune_plane_res_plot
  cmake.wall_generator
)

#set (FFLAGS "-qopenmp")
#set (LINK_FLAGS "-qopenmp")

if($ENV{ACC_ENABLE_MPI})
  list(APPEND EXE_SPECS "cmake.moga")
  list(APPEND EXE_SPECS "cmake.aperture_by_tracking")
  list(APPEND EXE_SPECS "cmake.touschek_background")
  list(APPEND EXE_SPECS "cmake.long_term_tracking_mpi")
  list(APPEND EXE_SPECS "cmake.tune_scan_mpi")
else()
  list(APPEND EXE_SPECS "cmake.long_term_tracking")
  list(APPEND EXE_SPECS "cmake.tune_scan")
endif()

## Profiler flags - Note: This is not compatible with the -fPIC flag, so -pg can't be set for SHARED build.
# set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# set (FFLAGS "-pg")  # For profiling
# set (CFLAGS "-pg")

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
