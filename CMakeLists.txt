
#set (LIBNAME util_programs)
cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(INC_DIRS 
  ../include
  ../cfortran/include
  ../cesr_utils/include
)

set(SRC_DIRS
)

set (EXE_SPECS 
  bmad_to_autocad.cmake
  bmad_to_csrtrack.cmake
  bmad_to_mad_and_xsif.cmake
  compare_tracking_methods_plot.cmake
  element_attributes.cmake
  f77_to_f90.cmake
  lat_to_wolski_orm_format.cmake
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





