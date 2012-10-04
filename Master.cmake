#-----------------------------------------------------------
# Master cmake lists file
# Implements the ACC build system.  Called by boilerplate
# text:
#      include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
# found in CMakeLists.txt files in project directories.
#-----------------------------------------------------------
cmake_minimum_required(VERSION 2.8)
#-----------------------------------------------------------
# Set link_directories relative path composition to use new
# behvaior that appends relative path information to the
# CMAKE_CURRENT_SOURCE_DIR value.
#-----------------------------------------------------------
cmake_policy(SET CMP0015 NEW)

#-------------------------------------------------------
# Import environment variables that influence the build
#-------------------------------------------------------
set(DISTRIBUTION_BUILD $ENV{DIST_BUILD})

IF (${DISTRIBUTION_BUILD})
  set (FORTRAN_COMPILER $ENV{DIST_F90})
ELSE ()
  set (FORTRAN_COMPILER "ifort")
ENDIF ()
  
IF (FORTRAN_COMPILER MATCHES "gfortran")

  set (RELEASE_NAME $ENV{DIST_BASE_DIR})
  set (RELEASE_NAME_TRUE "Off-site Distribution")
  set (RELEASE_DIR $ENV{DIST_BASE_DIR})
  set (COMPILER_CHOICE $ENV{DIST_F90})
  set (CMAKE_Fortran_COMPILER gfortran)
  set (COMPILER_SPECIFIC_F_FLAGS "-cpp -fno-range-check -fdollar-ok -fbacktrace -Bstatic -ffree-line-length-none")
  set (COMPILER_SPECIFIC_DEBUG_F_FLAGS "-O0")

ELSE ()

  set (RELEASE_NAME $ENV{ACC_RELEASE})
  set (RELEASE_NAME_TRUE $ENV{ACC_TRUE_RELEASE})
  set (RELEASE_DIR $ENV{ACC_RELEASE_DIR})
  set (CMAKE_Fortran_COMPILER ifort)
  set (COMPILER_SPECIFIC_F_FLAGS "-fpp")
  set (COMPILER_SPECIFIC DEBUG_F_FLAGS "-check bounds -check format -check uninit -warn declarations -ftrapuv")

ENDIF ()


set(PACKAGES_DIR ${RELEASE_DIR}/packages)

get_filename_component(SHORT_DIRNAME ${CMAKE_SOURCE_DIR} NAME)
#-------------------------------------------------------
# Check to see if shared object creation is enabled.
# If enabled, shared object libraries will only be
# created for projects that define 
#   CREATE_SHARED
# in their CmakeLists.txt file.
#
# TWO variables need to be set to produce a shared
# object library for a given project.
#  Shell: ACC_ENABLE_SHARED
#  Cmake: CREATE_SHARED
#-------------------------------------------------------
set(ENABLE_SHARED $ENV{ACC_ENABLE_SHARED})

IF (${LIBNAME})
	project(${LIBNAME})
ENDIF ()


#-----------------------------------
# C / C++ Compiler flags
#-----------------------------------
set (BASE_C_FLAGS "-Df2cFortran -O0 -std=gnu99 -mcmodel=medium -Wall -DCESR_UNIX -DCESR_LINUX -D_POSIX -D_REENTRANT -Wall -fPIC -Wno-trigraphs -Wno-unused")

set (BASE_CXX_FLAGS "-Wno-deprecated -mcmodel=medium -Wall -DCESR_UNIX -DCESR_LINUX -D_POSIX -D_REENTRANT -Wall -fPIC")


#-----------------------------------
# Fortran Compiler flags
#-----------------------------------
enable_language( Fortran )
set (BASE_Fortran_FLAGS "-Df2cFortran -DCESR_UNIX -DCESR_LINUX -u -traceback -mcmodel=medium ${COMPILER_SPECIFIC_F_FLAGS}")



set(CMAKE_EXE_LINKER_FLAGS "-lreadline -ltermcap -lcurses -lpthread -lstdc++ -lX11")


#--------------------------------------
# Honor requests for debug builds 
# composed with any variation in case.
#--------------------------------------
IF (CMAKE_BUILD_TYPE MATCHES "[Dd][Ee][Bb][Uu][Gg]")
  SET(DEBUG 1)
ELSE ()
  SET(DEBUG 0)
ENDIF ()


#------------------------------------------
# Honor requests for executable building
# made via environment variable.
#------------------------------------------
IF ($ENV{ACC_BUILD_EXES})
  SET(BUILD_EXES 1)
ELSE ()
  SET(BUILD_EXES 0)
ENDIF ()


#-----------------------------------------
# Print some friendly build information
# and according to the build type, set
# the following:
#    C Flags
#    F90 Flags
#-----------------------------------------
message("")
IF (DEBUG)
  message("Build type           : Debug")
  set (OUTPUT_BASEDIR ${CMAKE_SOURCE_DIR}/../debug)
  set (RELEASE_OUTPUT_BASEDIR ${RELEASE_DIR}/debug)
  set (PACKAGES_OUTPUT_BASEDIR ${RELEASE_DIR}/packages/debug)
  set(BASE_C_FLAGS "${BASE_C_FLAGS}")
  set(BASE_Fortran_FLAGS "${BASE_Fortran_FLAGS} ${COMPILER_SPECIFIC_DEBUG_F_FLAGS}")
ELSE ()
  message("Build type           : Production")
  set (OUTPUT_BASEDIR ${CMAKE_SOURCE_DIR}/../production)
  set (RELEASE_OUTPUT_BASEDIR ${RELEASE_DIR}/production)
  set (PACKAGES_OUTPUT_BASEDIR ${RELEASE_DIR}/packages/production)
  set (BASE_C_FLAGS "${BASE_C_FLAGS} -O2")
ENDIF ()
message("Linking with release : ${RELEASE_NAME} \(${RELEASE_NAME_TRUE}\)")
message("C Compiler           : ${CMAKE_C_COMPILER}")
message("Fortran Compiler     : ${CMAKE_Fortran_COMPILER}\n")


#-----------------------------------
# Output path definitions
#-----------------------------------
set (LIBRARY_OUTPUT_PATH ${OUTPUT_BASEDIR}/lib)
set (EXECUTABLE_OUTPUT_PATH ${OUTPUT_BASEDIR}/bin)
set (CMAKE_Fortran_MODULE_DIRECTORY ${OUTPUT_BASEDIR}/modules)


#-----------------------------------
#   System library locators
# (Also provides include directory
# locations.)
#-----------------------------------
find_package(X11)


#-------------------------
#   Include directories
#-------------------------
# TODO: Double each include directory entry to search for a local (../<xxx>) path and THEN
#       to a release-based path?
#
#  This leaves the possibility that someone may delete the local library, and initiate a build
#    that requires that library while leaving the local souce tree and include files intact.
#  This new build will then perform the divergent action of linking against the release library
# but extracting constants and other header information from the LOCAL source tree.  
#
SET (MASTER_INC_DIRS
  ${X11_INCLUDE_DIR}
  ${INC_DIRS}
  ${PACKAGES_DIR}/include
  ${PACKAGES_DIR}/forest/include
  ${PACKAGES_DIR}/recipes_c-ansi/include
  ${PACKAGES_DIR}/cfortran/include
  ${PACKAGES_DIR}/root/include
  ${PACKAGES_OUTPUT_BASEDIR}/modules
  ${ROOT_INC}
  ${RELEASE_DIR}/include
  ${OUTPUT_BASEDIR}/modules
  ${RELEASE_OUTPUT_BASEDIR}/modules
  ${RELEASE_DIR}/modules
)


#------------------------------------------------------
# Add local include paths to search list if they exist
#------------------------------------------------------
foreach(dir ${INC_DIRS})
  STRING(FIND ${dir} "../" relative)
  STRING(REPLACE "../" "" dirnew ${dir})
  IF (${relative} EQUAL 0)
    LIST(APPEND MASTER_INC_DIRS ${dir})
    LIST(APPEND MASTER_INC_DIRS ${RELEASE_DIR}/${dirnew})
  ELSE ()
    LIST(APPEND MASTER_INC_DIRS ${dir})
  ENDIF ()
endforeach(dir)


LIST(REMOVE_DUPLICATES MASTER_INC_DIRS)
INCLUDE_DIRECTORIES(${MASTER_INC_DIRS})


#-----------------------------------
#  Link directories - order matters
# Lowest level to highest, i.e. in
# order of increasing abstraction.
#-----------------------------------
SET(MASTER_LINK_DIRS
  /lib64
  /usr/lib64
  ${OUTPUT_BASEDIR}/lib
  ${PACKAGES_OUTPUT_BASEDIR}/lib
  ${PACKAGES_DIR}/root/lib
  ${RELEASE_OUTPUT_BASEDIR}/lib
  ${RELEASE_DIR}/lib
)
LINK_DIRECTORIES(${MASTER_LINK_DIRS})


#-------------------------------------------
# Collect list of all source files for all
# supported languages from all directories
# mentioned in project CMakeLists.txt file.
#-------------------------------------------
foreach(dir ${SRC_DIRS})

    set(temp_list)
    file(GLOB temp_list ${dir}/*.c)
    LIST(APPEND c_sources ${temp_list})

    file(GLOB temp_list ${dir}/*.C)
    LIST(APPEND c_sources ${temp_list})
    #-----------------------------------
    # Set compiler flag properties for 
    # all C source files.
    #-----------------------------------
    foreach (file ${c_sources})
      set_source_files_properties(${file} PROPERTIES COMPILE_FLAGS "${BASE_C_FLAGS} ${CFLAGS}")
    endforeach()
    LIST(APPEND sources ${c_sources})


    set(temp_list)
    file(GLOB temp_list ${dir}/*.cpp)
    LIST(APPEND cpp_sources ${temp_list})

    file(GLOB temp_list ${dir}/*.cc)
    LIST(APPEND cpp_sources ${temp_list})

    file(GLOB temp_list ${dir}/*.cxx)
    LIST(APPEND cpp_sources ${temp_list})
    #-----------------------------------
    # Set compiler flag properties for 
    # all C source files.
    #-----------------------------------
    foreach (file ${cpp_sources})
      set_source_files_properties(${file} PROPERTIES COMPILE_FLAGS "${BASE_CPP_FLAGS} ${CPPFLAGS}")
    endforeach()
    LIST(APPEND sources ${cpp_sources})


    set(temp_list)
    file(GLOB temp_list ${dir}/*.f)
    LIST(APPEND f_sources ${temp_list})

    file(GLOB temp_list ${dir}/*.F)
    LIST(APPEND f_sources ${temp_list})

    file(GLOB temp_list ${dir}/*.f90)
    LIST(APPEND f_sources ${temp_list})
    #-----------------------------------
    # Set compiler flag properties for
    # all Fortran source files.
    #-----------------------------------
    foreach (file ${f_sources})
      set_source_files_properties(${file} PROPERTIES COMPILE_FLAGS "${BASE_Fortran_FLAGS} ${FFLAGS}")
    endforeach()
    LIST(APPEND sources ${f_sources})

endforeach(dir)


set(DEPS)


#----------------------------------------------------------------
# If any pre-build script is specified, run it before building
# any code.  The pre-build script may generate code or header
# files.
#----------------------------------------------------------------
IF (PREBUILD_ACTION)
  message("Executing pre-build action...")
  EXECUTE_PROCESS (COMMAND ${PREBUILD_ACTION}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )
ENDIF ()

set(TARGETS)

IF (LIBNAME)
  add_library( ${LIBNAME} STATIC ${sources} )
  LIST(APPEND TARGETS ${LIBNAME})
  SET_TARGET_PROPERTIES(${LIBNAME} PROPERTIES OUTPUT_NAME ${LIBNAME})
ENDIF ()


#----------------------------------------------------------------
# Copy the contents of a config directory to 
# <DIR>/../config/${LIBNAME} if one exists.
#----------------------------------------------------------------
IF (IS_DIRECTORY "../config")
  message("Copying config directory contents to ${CMAKE_SOURCE_DIR}/../config/${LIBNAME}...")
  file (MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/../config")
  EXECUTE_PROCESS (COMMAND cp -rfu ${CMAKE_SOURCE_DIR}/config ${CMAKE_SOURCE_DIR}/../config/${SHORT_DIRNAME})
ENDIF ()


#----------------------------------------------------------------
# For selectively producing shared object libraries (.so files).
#
# Shell variable ACC_ENABLE_SHARED needs to be set to 
#  "Y" or "true" or "1"
#    - AND - 
# set (CREATE_SHARED true) needs to be present in the individual
#  project's CMakeLists.txt file.
#----------------------------------------------------------------
message("SHARED DEPS: ${SHARED_DEPS}")
IF (ENABLE_SHARED AND CREATE_SHARED)
  add_library (${LIBNAME}-shared SHARED ${sources})
  LIST(APPEND TARGETS ${LIBNAME}-shared)
  SET_TARGET_PROPERTIES (${LIBNAME}-shared PROPERTIES OUTPUT_NAME ${LIBNAME})
  TARGET_LINK_LIBRARIES (${LIBNAME}-shared ${SHARED_DEPS})
ENDIF ()



#---------------------------------------------------------------
# Process each EXE build description file mentioned in the
# project's CMakeLists.txt file.
#---------------------------------------------------------------
foreach(exespec ${EXE_SPECS})

  set(CFLAGS)
  set(FFLAGS)
  set(CPPFLAGS)

  set(SRC_DIRS)
  set(c_sources)
  set(cpp_sources)
  set(f_sources)

  include(${exespec})

  # TODO: Convert this to macro or function?
  foreach(dir ${INC_DIRS})
    STRING(FIND ${dir} "../" relative)
    STRING(REPLACE "../" "" dirnew ${dir})
    IF (${relative} EQUAL 0)
      LIST(APPEND MASTER_INC_DIRS ${dir})
      LIST(APPEND MASTER_INC_DIRS ${RELEASE_DIR}/${dirnew})
    ELSE ()
      LIST(APPEND MASTER_INC_DIRS ${dir})
    ENDIF ()
  endforeach(dir)
  LIST(REMOVE_DUPLICATES MASTER_INC_DIRS)
  INCLUDE_DIRECTORIES(${MASTER_INC_DIRS})

  set(DEPS ${LINK_LIBS})

  #----------------------------------------------------------------
  # Make this project's EXE build depend upon the product of each
  # build that is listed as a dependency.  If those binaries
  # are newer than this project's EXE, relink this project's EXE.
  # Only invoke add_library to tie in external dependencies a
  # single time for each unique target.
  #----------------------------------------------------------------
  foreach(dep ${DEPS})

    IF(${LIBNAME} MATCHES ${dep})
    ELSE()
      LIST(FIND TARGETS ${dep} DEP_SEEN)
      IF(${DEP_SEEN} EQUAL -1)
        IF (EXISTS ${OUTPUT_BASEDIR}/lib/lib${dep}.a)
          add_library(${dep} STATIC IMPORTED)
          LIST(APPEND TARGETS ${dep})
          set_property(TARGET ${dep} PROPERTY IMPORTED_LOCATION ${OUTPUT_BASEDIR}/lib/lib${dep}.a)
        ELSEIF (EXISTS ${PACKAGES_OUTPUT_BASEDIR}/lib/lib${dep}.a)
          add_library(${dep} STATIC IMPORTED)
          LIST(APPEND TARGETS ${dep})
          set_property(TARGET ${dep} PROPERTY IMPORTED_LOCATION ${PACKAGES_OUTPUT_BASEDIR}/lib/lib${dep}.a)
        ELSEIF (EXISTS ${RELEASE_DIR}/lib/lib${dep}.a)
          add_library(${dep} STATIC IMPORTED)
          LIST(APPEND TARGETS ${dep})
          set_property(TARGET ${dep} PROPERTY IMPORTED_LOCATION ${RELEASE_DIR}/lib/lib${dep}.a)
        ENDIF ()
      ENDIF()
    ENDIF()

  endforeach(dep)

  LINK_DIRECTORIES( ${LINK_DIRS} )

  #----------------------------------------------
  # Honor propagated control variable to build
  # any EXEs provided in a list from the main
  # CMakeLists.txt file.
  #----------------------------------------------
  IF (BUILD_EXES)
    set(BUILD_EXE_TOGGLE "")
  ELSE()
    set(BUILD_EXE_TOGGLE "EXCLUDE_FROM_ALL")
  ENDIF ()


  #-------------------------------------------
  # Collect list of all source files for all
  # supported languages from all directories
  # mentioned in project CMakeLists.txt file.
  #-------------------------------------------
  foreach(dir ${SRC_DIRS})

      set(temp_list)
      file(GLOB temp_list ${dir}/*.c)
      LIST(APPEND c_sources ${temp_list})

      file(GLOB temp_list ${dir}/*.C)
      LIST(APPEND c_sources ${temp_list})
      #-----------------------------------
      # Set compiler flag properties for 
      # all C source files.
      #-----------------------------------
      foreach (file ${c_sources})
	set_source_files_properties(${file} PROPERTIES COMPILE_FLAGS "${BASE_C_FLAGS} ${CFLAGS}")
      endforeach()
      LIST(APPEND SRC_FILES ${c_sources})


      set(temp_list)
      file(GLOB temp_list ${dir}/*.cpp)
      LIST(APPEND cpp_sources ${temp_list})

      file(GLOB temp_list ${dir}/*.cc)
      LIST(APPEND cpp_sources ${temp_list})

      file(GLOB temp_list ${dir}/*.cxx)
      LIST(APPEND cpp_sources ${temp_list})
      #-----------------------------------
      # Set compiler flag properties for 
      # all C source files.
      #-----------------------------------
      foreach (file ${cpp_sources})
	set_source_files_properties(${file} PROPERTIES COMPILE_FLAGS "${BASE_CPP_FLAGS} ${CPPFLAGS}")
      endforeach()
      LIST(APPEND SRC_FILES ${cpp_sources})


      set(temp_list)
      file(GLOB temp_list ${dir}/*.f)
      LIST(APPEND f_sources ${temp_list})

      file(GLOB temp_list ${dir}/*.F)
      LIST(APPEND f_sources ${temp_list})

      file(GLOB temp_list ${dir}/*.f90)
      LIST(APPEND f_sources ${temp_list})
      #-----------------------------------
      # Set compiler flag properties for
      # all Fortran source files.
      #-----------------------------------
      foreach (file ${f_sources})
	set_source_files_properties(${file} PROPERTIES COMPILE_FLAGS "${BASE_Fortran_FLAGS} ${FFLAGS}")
      endforeach()
      LIST(APPEND SRC_FILES ${f_sources})

  endforeach(dir)


  #----------------------------------------------
  # Apply user-specified compiler flags to each 
  # file being built into the executable.
  #----------------------------------------------
  SET(COMPILER_FLAGS "${BASE_C_FLAGS} ${COMPILER_FLAGS} ${CFLAGS}")
  SET(COMPILER_FLAGS "${BASE_Fortran_FLAGS} ${COMPILER_FLAGS} ${FFLAGS}")
  if (COMPILER_FLAGS)
    foreach(srcfile ${SRC_FILES})
      set(ext_match)
      STRING(FIND ${srcfile} ".c" ext_match)
      IF(NOT ${ext_match} EQUAL -1)
        set_source_files_properties(${srcfile} PROPERTIES COMPILE_FLAGS "${BASE_C_FLAGS} ${CFLAGS}")
      ENDIF()
      set(ext_match)
      STRING(FIND ${srcfile} ".C" ext_match)
      IF(NOT ${ext_match} EQUAL -1)
        set_source_files_properties(${srcfile} PROPERTIES COMPILE_FLAGS "${BASE_C_FLAGS} ${CFLAGS}")
      ENDIF()
      set(ext_match)
      STRING(FIND ${srcfile} ".cpp" ext_match)
      IF(NOT ${ext_match} EQUAL -1)
        set_source_files_properties(${srcfile} PROPERTIES COMPILE_FLAGS "${BASE_C_FLAGS} ${CFLAGS}")
      ENDIF()
      set(ext_match)
      STRING(FIND ${srcfile} ".cc" ext_match)
      IF(NOT ${ext_match} EQUAL -1)
        set_source_files_properties(${srcfile} PROPERTIES COMPILE_FLAGS "${BASE_C_FLAGS} ${CFLAGS}")
      ENDIF()
      set(ext_match)
      STRING(FIND ${srcfile} ".cxx" ext_match)
      IF(NOT ${ext_match} EQUAL -1)
        set_source_files_properties(${srcfile} PROPERTIES COMPILE_FLAGS "${BASE_C_FLAGS} ${CFLAGS}")
      ENDIF()
      #------------
      set(ext_match)
      STRING(FIND ${srcfile} ".f90" ext_match)
      IF(NOT ${ext_match} EQUAL -1)
        set_source_files_properties(${srcfile} PROPERTIES COMPILE_FLAGS "${BASE_Fortran_FLAGS} ${FFLAGS}")
      ENDIF()
      set(ext_match)
      STRING(FIND ${srcfile} ".f" ext_match)
      IF(NOT ${ext_match} EQUAL -1)
        set_source_files_properties(${srcfile} PROPERTIES COMPILE_FLAGS "${BASE_Fortran_FLAGS} ${FFLAGS}")
      ENDIF()
      set(ext_match)
      STRING(FIND ${srcfile} ".F" ext_match)
      IF(NOT ${ext_match} EQUAL -1)
        set_source_files_properties(${srcfile} PROPERTIES COMPILE_FLAGS "${BASE_Fortran_FLAGS} ${FFLAGS}")
      ENDIF()
      set(ext_match)
    endforeach(srcfile)
  endif ()


  # Actually request creation of executable target
  add_executable(${EXENAME}-exe
      ${BUILD_EXE_TOGGLE}
      ${SRC_FILES}
  )

  SET_TARGET_PROPERTIES(${EXENAME}-exe
          PROPERTIES
          OUTPUT_NAME
          ${EXENAME}
  )


  # Create map file output directory if it doesn't yet exist.
  IF(IS_DIRECTORY ${OUTPUT_BASEDIR}/map)
  ELSE()
    file(MAKE_DIRECTORY ${OUTPUT_BASEDIR}/map)
  ENDIF()

  # Set up linking for the executable.
  # Always produce a map file.  It is placed in the ../<build_type>/map directory
  # created during build setup.
  TARGET_LINK_LIBRARIES(${EXENAME}-exe
          ${LINK_LIBS}
          ${LINK_FLAGS} "-Wl,-Map=${OUTPUT_BASEDIR}/map/${EXENAME}.map"
  )

  SET(CFLAGS)
  SET(FFLAGS)
  SET(COMPILER_FLAGS)
  SET(LINK_FLAGS)

endforeach(exespec)


#-------------------------------------------------------------------
# If a shared object library has been built, and if a Makefile.mex
# file exists for building Matlab MEX wrappers, call that makefile.
#-------------------------------------------------------------------
FOREACH(target ${TARGETS})

  IF(target MATCHES ${LIBNAME}-shared)
    IF(EXISTS ../Makefile.mex)
      ADD_CUSTOM_COMMAND(TARGET ${LIBNAME}-shared
      POST_BUILD
      COMMAND gmake -f Makefile.mex
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      )
    ENDIF(EXISTS ../Makefile.mex)
  ENDIF()

ENDFOREACH()
