cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

SET (LIBNAME tao)

SET (CREATE_SHARED true)

SET (PLOT_LINK_LIBS $ENV{PLOT_LINK_LIBS})

SET (SRC_DIRS
  code 
  ## hook # Now using function pointers
  version
  )

SET (EXE_SPECS
  cmake.tao
  )

SET (SHARED_DEPS
  bmad
  sim_utils
  ${PLOT_LINK_LIBS}
  $ENV{ACC_BMAD_LINK_LIBS}
  gomp
  X11
  readline
)

SET (TAO_SOURCE_ROOT_CMAKE "${CMAKE_CURRENT_SOURCE_DIR}/..")
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/code/tao_preprocess_path.f90.in
  ${CMAKE_CURRENT_BINARY_DIR}/tao_preprocess_path.f90
  @ONLY
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)

# Add the generated module source to the tao library target.
target_sources(tao PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/tao_preprocess_path.f90)
