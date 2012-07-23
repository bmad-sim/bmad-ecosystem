set (LIBNAME bsim)
cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(INC_DIRS 
  ../include
  synrad3d
)

#  /usr/lib64/openmpi/1.4-gcc/include
#  /usr/include/compat-openmpi-x86_64

set (SRC_DIRS
  code
)

set (DEPS
  bmadz
  cesr_utils
  bmad
  sim_utils
  recipes_f-90_LEPP
  forest
  pgplot
  xsif
)

set (EXE_SPECS
  analyzer.cmake
  bbu.cmake
  dynamic_aperture.cmake
  synrad.cmake
  synrad3d.cmake
  synrad_aperture_to_wall.cmake
  tune_plane_res_plot.cmake
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
