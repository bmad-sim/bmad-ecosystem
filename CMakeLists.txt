set (LIBNAME bsim)
cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

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
  #cmake.envelope_ibs
  cmake.srdt_lsq_soln
  cmake.long_term_tracking
  cmake.frequency_map
  #cmake.frequency_map_mpi
  cmake.bbu
  cmake.dynamic_aperture
  cmake.ibs_linac
  cmake.ibs_ring
  cmake.spin_stroboscope
  cmake.synrad
  cmake.synrad3d
  cmake.synrad_aperture_to_wall
  cmake.tune_plane_res_plot
)

#set (FFLAGS "-qopenmp")
#set (LINK_FLAGS "-qopenmp")

if($ENV{ACC_ENABLE_MPI})
  list(APPEND EXE_SPECS "cmake.moga")
  list(APPEND EXE_SPECS "cmake.long_term_tracking_mpi")
endif()

#set (FFLAGS "-pg")  # For profiling
#set (CFLAGS "-pg")

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
