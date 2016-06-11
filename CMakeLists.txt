set (LIBNAME util_programs)
cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(INC_DIRS 
  ../include
)

set(SRC_DIRS 
  modules
)

set (FFLAGS "-coarray")
set (LINK_FLAGS "-coarray")

set (EXE_SPECS 
  # cmake.moga
  cmake.wake_fit
  cmake.bmad_to_autocad
  cmake.bmad_to_csrtrack
  cmake.bmad_to_mad_and_xsif
  cmake.compare_tracking_methods_plot
  cmake.compare_tracking_methods_text
  cmake.element_attributes
  cmake.f77_to_f90
  cmake.lattice_cleaner
  cmake.sad_to_bmad_postprocess
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





