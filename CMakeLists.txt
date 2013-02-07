set (LIBNAME bsim)
cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(INC_DIRS 
  ../include
)

#  synrad3d
set (SRC_DIRS
  code
  synrad3d
)

set (EXE_SPECS
  analyzer.cmake
  bbu.cmake
  dynamic_aperture.cmake
  synrad.cmake
  synrad3d.cmake
  synrad_aperture_to_wall.cmake
  tune_plane_res_plot.cmake
  find_quads.cmake
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
