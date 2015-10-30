set (LIBNAME bsim)
cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(INC_DIRS 
)

set (SRC_DIRS
  code
  synrad_lib
)

set (EXE_SPECS
  cmake.s3drefl
  cmake.bbu
  cmake.bbu_old
  cmake.dynamic_aperture
  cmake.ibs_linac
  cmake.ibs_ring
  cmake.synrad
  cmake.synrad3d
  cmake.synrad_aperture_to_wall
  cmake.tune_plane_res_plot
)

#set (FFLAGS "-openmp")
#set (LINK_FLAGS "-openmp")

#set (FFLAGS "-pg")  # For profiling
#set (CFLAGS "-pg")

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
