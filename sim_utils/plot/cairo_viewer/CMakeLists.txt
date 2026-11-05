# qp_cairo_viewer: the on-screen window for the quick_plot Cairo backend.
#
# Built automatically as part of sim_utils when ACC_PLOT_PACKAGE = "cairo".
# Can also be built standalone:
#   cmake -S . -B build && cmake --build build
#
# Needs SDL2 and Cairo (MacPorts: libsdl2 cairo; conda-forge: sdl2 cairo; MSYS2: SDL2 cairo).
# If SDL2 is not found the viewer is skipped with a warning: hardcopy output still works
# but on-screen ("X") plot pages cannot be opened.

cmake_minimum_required(VERSION 3.13)
project(qp_cairo_viewer C)

find_package(SDL2 QUIET)
if (NOT SDL2_FOUND)
  message(WARNING "SDL2 not found: qp_cairo_viewer (the quick_plot Cairo plot window) will not be built.")
  return()
endif ()

find_package(PkgConfig REQUIRED)
pkg_check_modules(QP_VIEWER_CAIRO REQUIRED IMPORTED_TARGET cairo)

add_executable(qp_cairo_viewer qp_cairo_viewer.c ../qp_cairo_render.c)
target_compile_definitions(qp_cairo_viewer PRIVATE QP_CAIRO_VIEWER)
target_include_directories(qp_cairo_viewer PRIVATE ..)

if (TARGET SDL2::SDL2)
  target_link_libraries(qp_cairo_viewer PRIVATE SDL2::SDL2 PkgConfig::QP_VIEWER_CAIRO)
else ()
  target_include_directories(qp_cairo_viewer PRIVATE ${SDL2_INCLUDE_DIRS})
  target_link_libraries(qp_cairo_viewer PRIVATE ${SDL2_LIBRARIES} PkgConfig::QP_VIEWER_CAIRO)
endif ()

if (NOT WIN32)
  target_link_libraries(qp_cairo_viewer PRIVATE m)
endif ()

install(TARGETS qp_cairo_viewer RUNTIME DESTINATION bin)
