if (DEFINED ENV{BUILD_TYPE})
  cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
  project(ACC)

  SET (LIBNAME tao)

  SET (CREATE_SHARED true)

  SET (PLOT_LINK_LIBS $ENV{PLOT_LINK_LIBS})

  SET (SRC_DIRS
    code
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

  include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
else()
  file(GLOB src code/*.f90)

  add_library(tao SHARED ${src})
  target_link_libraries(tao PUBLIC bmad)
  set_target_properties(tao PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
  )
  set_source_files_properties(
    code/tao_get_user_input_mod.f90
    PROPERTIES Fortran_PREPROCESS ON)
  target_include_directories(tao PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
    $<INSTALL_INTERFACE:modules>
  )
  target_link_libraries(tao PUBLIC bmad_options bmad sim_utils)

  add_executable(tao_exe program/tao_program.f90)
  set_target_properties(tao_exe PROPERTIES OUTPUT_NAME tao)
  target_link_libraries(tao_exe tao)

  install(TARGETS tao EXPORT Bmad)
  install(TARGETS tao_exe)
  foreach(d code program)
    install(DIRECTORY ${d} DESTINATION src/tao
      FILES_MATCHING
      PATTERN "*.f"
      PATTERN "*.f90"
      PATTERN "*.c"
      PATTERN "*.cc"
      PATTERN "*.cpp")
  endforeach()
  install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/" DESTINATION modules FILES_MATCHING PATTERN "*.mod")
endif()
