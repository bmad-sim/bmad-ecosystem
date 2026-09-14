cmake_minimum_required(VERSION 3.14)
project(Bmad LANGUAGES Fortran C CXX)

option(ENABLE_OPENMP "Enable OpenMP support" OFF)
if (ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True)
set(CMAKE_BUILD_WITH_INSTALL_PATH True)
if(APPLE)
  set(CMAKE_MACOSX_RPATH ON)
  set(CMAKE_INSTALL_RPATH "@loader_path../lib")
else()
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
endif()

find_package(PkgConfig REQUIRED)

add_library(bmad_options INTERFACE)
target_compile_options(bmad_options INTERFACE
  $<$<COMPILE_LANGUAGE:Fortran>:-fdollar-ok -ffree-line-length-none>)
install(TARGETS bmad_options EXPORT Bmad)
add_subdirectory(forest)
add_subdirectory(sim_utils)
add_subdirectory(bmad)
add_subdirectory(tao)
add_subdirectory(cpp_bmad_interface)
add_subdirectory(bsim)
install(EXPORT Bmad NAMESPACE Bmad:: DESTINATION lib/cmake/Bmad)
include(CMakePackageConfigHelpers)
configure_package_config_file(BmadConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/BmadConfig.cmake
  INSTALL_DESTINATION lib/cmake/Bmad)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/BmadConfig.cmake DESTINATION lib/cmake/Bmad)
