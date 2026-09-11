if (DEFINED ENV{BUILD_TYPE})
  cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
  project(ACC)

  set (LIBNAME cpp_bmad_interface)

  set(INC_DIRS 
    include
  )

  set(SRC_DIRS
    code
  )

  set(EXE_SPECS 
    cmake.cpp_bmad_interface_test
  )

  include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
else()
  add_library(cpp_bmad_interface SHARED code/bmad_cpp_convert_mod.f90 code/cpp_bmad_convert.cpp code/cpp_equality.cpp)
  set_target_properties(cpp_bmad_interface PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
  )
  target_link_libraries(cpp_bmad_interface PUBLIC bmad_options bmad)
  target_include_directories(cpp_bmad_interface PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:modules>
  )

  install(TARGETS cpp_bmad_interface EXPORT Bmad)
  install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
  install(DIRECTORY code DESTINATION src/cpp_bmad_interface FILES_MATCHING PATTERN "*.cpp" PATTERN "*.f90")
  install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/" DESTINATION modules FILES_MATCHING PATTERN "*.mod")
endif()
