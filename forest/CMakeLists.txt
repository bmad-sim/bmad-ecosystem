if (DEFINED ENV{ACC_SET_GMAKE_JOBS} AND DEFINED ENV{GMAKE})
  cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
  project(ACC)

  SET (LIBNAME forest)

  SET (CREATE_SHARED true)

  SET (INC_DIRS 
    include
  )

  SET (SRC_DIRS
    code
  )

  include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
else()
  file(GLOB forest_src code/*.f90)
  add_library(forest SHARED ${forest_src})

  set_target_properties(forest PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
  )

  target_link_libraries(forest PUBLIC bmad_options)

  target_include_directories(forest PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:modules>
  )

  install(TARGETS forest EXPORT Bmad)
  install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.inc")
  install(DIRECTORY code DESTINATION src/forest FILES_MATCHING PATTERN "*.f90")
  install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/" DESTINATION modules FILES_MATCHING PATTERN "*.mod")
endif()
