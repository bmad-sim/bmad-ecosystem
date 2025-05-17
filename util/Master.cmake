#-----------------------------------------------------------
# Master CMake file.
# Implements the ACC build system.
# How to include:
#      include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
# Found in CMakeLists.txt files in project directories.

#-----------------------------------------------------------
# Set link_directories relative path composition to use new
# behvaior that appends relative path information to the
# CMAKE_CURRENT_SOURCE_DIR value.
#-----------------------------------------------------------
cmake_policy (SET CMP0015 NEW)


#-----------------------------------------------------------
# Set CESR_FLAGS depening on OS type
#-----------------------------------------------------------

IF (${WIN32})
    SET (CESR_FLAGS "-DCESR_WINCVF")
ELSE ()
    SET (CESR_FLAGS "-DCESR_UNIX -DCESR_LINUX")
ENDIF ()

#------------------------------------------
# Honor requests for compiling with openmp
# made via environment variable.
#------------------------------------------
IF ($ENV{ACC_ENABLE_OPENMP})
  SET (ACC_ENABLE_OPENMP 1)
ELSE ()
  SET (ACC_ENABLE_OPENMP 0)
ENDIF ()


#------------------------------------------
# Honor requests for compiling with OpenMPI
# made via environment variable.
#------------------------------------------
IF ($ENV{ACC_ENABLE_MPI})
  SET (ACC_ENABLE_MPI 1)
ELSE ()
  SET (ACC_ENABLE_MPI 0)
ENDIF ()


#------------------------------------------
# Honor requests for gfortran compiling with 
# -O2 flag made via environment variable.
#------------------------------------------
IF ($ENV{ACC_ENABLE_GFORTRAN_OPTIMIZATION}) 
  IF ("$ENV{DIST_F90}" MATCHES "gfortran")
    SET (ACC_GFORTRAN_OPTIMIZATION_FLAG "-O2")
  ELSEIF ("$ENV{ACC_SET_F_COMPILER}" MATCHES "gfortran")
    SET (ACC_GFORTRAN_OPTIMIZATION_FLAG "-O2")
  ELSE ()
    SET (ACC_GFORTRAN_OPTIMIZATION_FLAG)
  ENDIF ()
ENDIF ()


#-------------------------------------------------------
# Import environment variables that influence the build
# Also disable the inclusion of RPATH for Releases but not Distributions RT#64118
#-------------------------------------------------------

set (FORTRAN_COMPILER $ENV{DIST_F90})
set (RELEASE_DIR $ENV{DIST_BASE_DIR})
IF (EXISTS ${RELEASE_DIR}/packages)
  set (PACKAGES_DIR ${RELEASE_DIR}/packages)
ELSE ()
  set (PACKAGES_DIR ${RELEASE_DIR})
ENDIF ()

# Explicitly remove 32-bit Libraries from Linux build PATH for 64-bit builds - RT#43178
# Added /lib to further restict 32-bit Linux build PATH - RT#56203
IF (${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND NOT "$ENV{ACC_FORCE_32_BIT}" MATCHES "Y")
  SET (CMAKE_IGNORE_PATH /lib /usr/lib)
  SET (CMAKE_SYSTEM_IGNORE_PATH /lib /usr/lib)
ENDIF ()

  
IF (FORTRAN_COMPILER MATCHES "gfortran")

  SET (RELEASE_NAME $ENV{DIST_BASE_DIR})
  IF ("$ENV{DIST_F90}" MATCHES "gfortran")
    SET (RELEASE_NAME_TRUE "Off-site Distribution")
  ENDIF ()
  SET (COMPILER_CHOICE ${FORTRAN_COMPILER})
  SET (CMAKE_Fortran_COMPILER gfortran)
     IF (${ACC_ENABLE_OPENMP})
       SET (COMPILER_SPECIFIC_F_FLAGS "-cpp -fno-range-check -fdollar-ok -fbacktrace -Bstatic -ffree-line-length-none -fopenmp")
       SET (OPENMP_LINK_LIBS "gomp")
       IF($ENV{ACC_ENABLE_SHARED} OR $ENV{ACC_ENABLE_SHARED_ONLY})
         SET (SHARED_DEPS ${OPENMP_LINK_LIBS} ${SHARED_DEPS})
         SET (CMAKE_SHARED_LINKER_FLAGS "-fopenmp ${CMAKE_SHARED_LINKER_FLAGS}")
       ENDIF ()
     ELSE ()
       SET (COMPILER_SPECIFIC_F_FLAGS "-cpp -fno-range-check -fdollar-ok -fbacktrace -Bstatic -ffree-line-length-none")
     ENDIF ()
  SET (COMPILER_SPECIFIC_DEBUG_F_FLAGS "-O0 -fno-range-check -fbounds-check -Wuninitialized")

ELSE ()

  SET (RELEASE_NAME $ENV{ACC_RELEASE})
  SET (RELEASE_NAME_TRUE $ENV{ACC_TRUE_RELEASE})
  SET (CMAKE_Fortran_COMPILER ifort)
     IF (${ACC_ENABLE_OPENMP})
       SET (OPENMP_LINK_LIBS "-liomp5")
       EXECUTE_PROCESS (
         COMMAND bash -c "ifort --version | head -1 | awk ' { print $3 } '"
         OUTPUT_VARIABLE INTEL_VERSION_OUPUT
         )
       IF ( ${INTEL_VERSION_OUPUT} VERSION_GREATER "16.0.0" )
         SET (IFORT_OPENMP_FLAG "-qopenmp")
       ELSE()
         SET (IFORT_OPENMP_FLAG "-openmp")
       ENDIF ()
       SET (COMPILER_SPECIFIC_F_FLAGS "-fpp ${IFORT_OPENMP_FLAG}")
       IF($ENV{ACC_ENABLE_SHARED} OR $ENV{ACC_ENABLE_SHARED_ONLY})
         SET (SHARED_DEPS ${OPENMP_LINK_LIBS} ${SHARED_DEPS})
       ENDIF ()
     ELSE ()
       SET (COMPILER_SPECIFIC_F_FLAGS "-fpp")
     ENDIF ()
  SET (COMPILER_SPECIFIC_DEBUG_F_FLAGS "-check bounds -check format -check uninit -warn declarations -ftrapuv")

ENDIF ()

IF (${ACC_ENABLE_MPI})
  SET (CMAKE_Fortran_COMPILER mpifort)
  EXEC_PROGRAM (mpifort ARGS --showme:compile OUTPUT_VARIABLE MPI_COMPILE_FLAGS)
  EXEC_PROGRAM (mpifort ARGS --showme:link OUTPUT_VARIABLE MPI_LINK_FLAGS)
  EXEC_PROGRAM (mpifort ARGS --showme:incdirs OUTPUT_VARIABLE MPI_INC_DIR)
  EXEC_PROGRAM (mpifort ARGS --showme:libdirs OUTPUT_VARIABLE MPI_LIB_DIR)
  EXEC_PROGRAM (mpifort ARGS --showme:libs OUTPUT_VARIABLE MPI_LIBS)
  EXEC_PROGRAM (mpifort ARGS --showme:version OUTPUT_VARIABLE MPI_FORT_VERSION)
  # The preprocessor option flag "-DACC_MPI" was requested and tested in RT#31494 but not implemented until RT#58944
  # The -fallow-argument-mismatch flag allows, for mpi routines, actual argument type mismatch which is not an error.
  IF (${CMAKE_Fortran_COMPILER} STREQUAL "gfortran")
    SET (MPI_COMPILE_FLAGS "-DACC_MPI -fallow-argument-mismatch ${MPI_COMPILE_FLAGS}")
  ELSE ()
    SET (MPI_COMPILE_FLAGS "-DACC_MPI ${MPI_COMPILE_FLAGS}")
  ENDIF ()
ENDIF ()

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
SET (ENABLE_SHARED $ENV{ACC_ENABLE_SHARED})

# Added support for SHARED only library builds, as requested in RT#63875
SET (ENABLE_SHARED_ONLY $ENV{ACC_ENABLE_SHARED_ONLY})
IF (${ENABLE_SHARED_ONLY}) # Force setting ENABLE_SHARED
  SET (ENABLE_SHARED "Y")
ENDIF ()

IF (${LIBNAME})
  PROJECT (${LIBNAME})
ENDIF ()


#-----------------------------------
#   System library locators
# (Also provides include directory
# locations.)
#-----------------------------------

find_package(X11)

#-----------------------------------
# Plotting Library Linker flags
#-----------------------------------

SET (PLOT_LINK_LIBS $ENV{PLOT_LINK_LIBS})

IF ($ENV{ACC_PLOT_PACKAGE} MATCHES "plplot")
  SET (PLOT_LIBRARY_FLAG "-DCESR_PLPLOT")

  IF (NOT "$ENV{ACC_PLOT_DISPLAY_TYPE}" MATCHES "QT")
    SET (PLOT_LINK_FLAGS "-lX11 $ENV{PLOT_LINK_FLAGS}")
  ELSE ()
    SET (PLOT_LINK_FLAGS "$ENV{PLOT_LINK_FLAGS}")
  ENDIF ()

  IF (${MSYS})
    # Assuming Qt backend:
    SET (SHARED_LINK_LIBS QtSvg4 QtGui4 QtCore4 ${SHARED_LINK_LIBS})
  ELSE ()
    SET (SHARED_LINK_LIBS cairo pango-1.0 pangocairo-1.0 gobject-2.0 ${SHARED_LINK_LIBS})
  ENDIF ()

  SET (ACC_PLOT_INC_DIRS)

ELSEIF ($ENV{ACC_PLOT_PACKAGE} MATCHES "pgplot")
  SET (PLOT_LIBRARY_FLAG "-DCESR_PGPLOT")
  IF ("$ENV{ACC_ENABLE_SHARED}" MATCHES "Y")
# This is a temporary fix -amd275 
#    SET (PLOT_LINK_FLAGS "-lX11 $ENV{PLOT_LINK_FLAGS}")
    SET (PLOT_LINK_FLAGS "-lX11")
  ELSE ()
    SET (PLOT_LINK_FLAGS "-lX11")
  ENDIF ()

ELSEIF ($ENV{ACC_PLOT_PACKAGE} MATCHES "none")
  SET (PLOT_LIBRARY_FLAG "-DCESR_NOPLOT")
  SET (PLOT_LINK_LIBS)
  SET (PLOT_LINK_FLAGS)

ENDIF ()

IF (${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND NOT $ENV{ACC_PLOT_PACKAGE} MATCHES "none")
  SET (ACC_PLOT_LIB_DIRS /usr/lib64)
ENDIF ()

IF (${CMAKE_Fortran_COMPILER} STREQUAL "ifort" AND "$ENV{ACC_ENABLE_SHARED}" MATCHES "Y")
  SET (PLOT_LINK_FLAGS "${PLOT_LINK_FLAGS} -lifcore -lifport -limf -lsvml -lintlc")
ENDIF ()

#-----------------------------------
# C / C++ Compiler flags
#-----------------------------------

SET (BASE_C_FLAGS "-Df2cFortran -O0 -std=gnu99 ${CESR_FLAGS} ${PLOT_LIBRARY_FLAG} -D_POSIX -D_REENTRANT -Wall -fPIC -Wno-trigraphs -Wno-unused ${MPI_COMPILE_FLAGS}")
SET (BASE_CXX_FLAGS "-O0 -Wno-deprecated ${CESR_FLAGS} ${PLOT_LIBRARY_FLAG} -D_POSIX -D_REENTRANT -Wall -fPIC -Wno-trigraphs -Wno-unused ${MPI_COMPILE_FLAGS}")

#-----------------------------------                                                                                
# For non-Linux or non-ifort or 
# non-64-bit builds, do not use
# unspported flag option 
# "-mcmodel=medium"                                                                      
#----------------------------------- 
                                                                               
IF (${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND ${FORTRAN_COMPILER} MATCHES "ifort" AND CMAKE_SIZEOF_VOID_P EQUAL 8)
  SET (BASE_C_FLAGS "${BASE_C_FLAGS} -mcmodel=medium")
  SET (BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -mcmodel=medium")
ENDIF ()

#--------------------------------------
# Set STDCXX_LINK_LIBS and STDCXX_LINK_FLAGS variables
# As defined in build_flags_config
# See RT#58127
#--------------------------------------

SET (STDCXX_LINK_LIBS $ENV{STDCXX_LINK_LIBS})
SET (STDCXX_LINK_FLAGS $ENV{STDCXX_LINK_FLAGS})

#-----------------------------------
# Readline Library Definitions 
#-----------------------------------

SET (READLINE_LINK_LIBS readline)
SET (READLINE_LINK_FLAGS "-lreadline")

#--------------------------------------
# Fortran Compiler flags
#--------------------------------------

enable_language( Fortran )

SET (BASE_Fortran_FLAGS "-Df2cFortran ${CESR_FLAGS} -u -traceback ${COMPILER_SPECIFIC_F_FLAGS} ${PLOT_LIBRARY_FLAG} ${MPI_COMPILE_FLAGS}")

# For ifort, move STACK memory to the HEAP, see RT#45710
IF (${FORTRAN_COMPILER} MATCHES "ifort")
  SET (BASE_Fortran_FLAGS "${BASE_Fortran_FLAGS} -heap-arrays 32")
ENDIF ()

IF (${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND ${FORTRAN_COMPILER} MATCHES "ifort" AND CMAKE_SIZEOF_VOID_P EQUAL 8)
  SET (BASE_Fortran_FLAGS "${BASE_Fortran_FLAGS} -mcmodel=medium")
ENDIF ()

IF ($ENV{ACC_ENABLE_FPIC})
   SET (BASE_Fortran_FLAGS "${BASE_Fortran_FLAGS} -fPIC")
ENDIF ()

IF ($ENV{ACC_ENABLE_SHARED})
   SET (BASE_Fortran_FLAGS "${BASE_Fortran_FLAGS} -fPIC")
ENDIF ()

    SET (ACC_LINK_FLAGS ${ACC_LINK_FLAGS} ${MPI_LINK_FLAGS} ${PLOT_LINK_FLAGS} ${STDCXX_LINK_FLAGS})

if (FORTRAN_COMPILER MATCHES "gfortran")
 list (APPEND ACC_LINK_FLAGS "-ldl")
endif()

IF (${MSYS})
    SET (ACC_LINK_FLAGS)
    SET (SHARED_LINK_LIBS ${SHARED_LINK_LIBS} ${STDCXX_LINK_LIBS} readline termcap gdi32 Comdlg32)
ENDIF ()

string (STRIP "${ACC_LINK_FLAGS}" ACC_LINK_FLAGS)

IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND NOT $ENV{ACC_CONDA_BUILD} MATCHES "Y")
  find_program(BREW "brew")
  find_program(PORT "port")
  set(DO_BREW "NO")
  set(DO_PORT "NO")
  if (NOT BREW STREQUAL "BREW-NOTFOUND")
    set(DO_BREW "YES")
  endif()
  if (NOT PORT STREQUAL "PORT-NOTFOUND")
    set(DO_PORT "YES")
  endif()
  if (DO_BREW AND DO_PORT)
    set(DO_BREW "NO")
    set(DO_PORT "NO")
    if (DEFINED ENV{BMAD_MAC_PACKAGE})
      string(TOLOWER "$ENV{BMAD_MAC_PACKAGE}" MAC_PACKAGE)
      if (MAC_PACKAGE STREQUAL "macports")
	set (DO_PORT "YES")
      elseif (MAC_PACKAGE STREQUAL "homebrew")
	set (DO_BREW "YES")
      endif ()
    endif()
    if (NOT (DO_BREW OR DO_PORT))
      message(FATAL_ERROR
	"Both HomeBrew and Macports are installed.\n"
	"Set BMAD_MAC_PACKAGE to either \"macports\" or \"homebrew\" to choose which to use")
    endif()
  endif()
  if (DO_BREW)
    execute_process(COMMAND "brew" "--prefix" OUTPUT_VARIABLE HOMEBREW_PREFIX)
    string(STRIP "${HOMEBREW_PREFIX}" HOMEBREW_PREFIX)
    set(ACC_INC_DIRS ${ACC_INC_DIRS} "${HOMEBREW_PREFIX}/include" "${HOMEBREW_PREFIX}/opt/readline/include")
    set(ACC_LIB_DIRS ${ACC_LIB_DIRS} "${HOMEBREW_PREFIX}/lib" "${HOMEBREW_PREFIX}/opt/readline/lib")
  elseif (DO_PORT)
    set(ACC_INC_DIRS ${ACC_INC_DIRS} /opt/local/include)
    set(ACC_LIB_DIRS ${ACC_LIB_DIRS} /opt/local/lib)
  else()
    message(FATAL_ERROR "This is MacOS, but neither HomeBrew nor MacPorts is installed")
  endif()
ENDIF()

# Environment variables BMAD_USER_INC_DIRS and BMAD_USER_LIB_DIRS hold include/library directories
# to search for include/module or libraries respectively. Multiple directories separated by semicolons.
SET (ACC_INC_DIRS $ENV{BMAD_USER_INC_DIRS} ${ACC_PLOT_INC_DIRS} ${ACC_INC_DIRS} ${MPI_INC_DIRS})
SET (ACC_LIB_DIRS $ENV{BMAD_USER_LIB_DIRS} ${ACC_PLOT_LIB_DIRS} ${ACC_LIB_DIRS} ${MPI_LIB_DIRS})


#--------------------------------------
# Set ACC_BMAD_LINK_LIBS variable
# As requested by dcs16 in RT#56661
#--------------------------------------

SET (ACC_BMAD_LINK_LIBS $ENV{ACC_BMAD_LINK_LIBS})


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


#------------------------------------------
# Honor requests for test executable 
# building made via environment variable.
#------------------------------------------
IF ($ENV{ACC_BUILD_TEST_EXES})
  SET(BUILD_TEST_EXES 1)
ELSE ()
  SET(BUILD_TEST_EXES 0)
ENDIF ()


#-----------------------------------------
# Print some friendly build information
# and according to the build type, set
# the following:
#    C Flags
#    F90 Flags
#-----------------------------------------
message("")
execute_process(COMMAND date OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE this_date)

IF (DEBUG)
  message("Build type           : Debug   -- ${this_date}")
  set (OUTPUT_BASEDIR ${CMAKE_SOURCE_DIR}/../debug)
  set (RELEASE_OUTPUT_BASEDIR ${RELEASE_DIR}/debug)
  set (PACKAGES_OUTPUT_BASEDIR ${PACKAGES_DIR}/debug)
  set (BASE_Fortran_FLAGS "${BASE_Fortran_FLAGS} ${COMPILER_SPECIFIC_DEBUG_F_FLAGS}")
ELSE ()
  message("Build type           : Production   -- ${this_date}")
  set (OUTPUT_BASEDIR ${CMAKE_SOURCE_DIR}/../production)
  set (RELEASE_OUTPUT_BASEDIR ${RELEASE_DIR}/production)
  set (PACKAGES_OUTPUT_BASEDIR ${PACKAGES_DIR}/production)
  set (BASE_C_FLAGS "${BASE_C_FLAGS} -O2")
  set (BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -O2")
  set (BASE_Fortran_FLAGS "${BASE_Fortran_FLAGS} ${ACC_GFORTRAN_OPTIMIZATION_FLAG}")
ENDIF ()

message("Current Directory    : ${CMAKE_CURRENT_SOURCE_DIR}")
message("Linking with release : ${RELEASE_NAME} \(${RELEASE_NAME_TRUE}\)")
message("C Compiler           : ${CMAKE_C_COMPILER}")
message("Fortran Compiler     : ${CMAKE_Fortran_COMPILER}")
message("Plotting Libraries   : ${PLOT_LINK_LIBS}")

IF (DEFINED SHARED_LINK_LIBS)
message("Shared Libraries     : ${SHARED_LINK_LIBS}")
ENDIF()
IF ($ENV{ACC_ENABLE_OPENMP})
  IF (${FORTRAN_COMPILER} MATCHES "ifort")
    MESSAGE("OpenMP ifort Flag    : ${IFORT_OPENMP_FLAG}")
  ELSE ()
    message("OpenMP gfortran Flag : -fopenmp")
    message("OpenMP Linker Libs   : ${OPENMP_LINK_LIBS}")
  ENDIF()
ELSE ()
  message("OpenMP Support       : Not Enabled")
ENDIF ()
IF ($ENV{ACC_ENABLE_MPI})
  MESSAGE("MPI Support          : ${MPI_FORT_VERSION}")
  MESSAGE("MPI Compiler Flags   : ${MPI_COMPILE_FLAGS}")
  MESSAGE("MPI Link Flags       : ${MPI_LINK_FLAGS}")
ELSE()
  MESSAGE("MPI Support          : Not Enabled")
ENDIF()
message("FFLAGS               : ${FFLAGS} ${FCFLAGS}")
message("${FORTRAN_COMPILER} Compiler Flags : ${BASE_Fortran_FLAGS}")
message("${FORTRAN_COMPILER} Linker Flags   : ${ACC_LINK_FLAGS} ${OPENMP_LINK_LIBS}")

MESSAGE("ACC_BMAD_LINK_LIBS Link Library list: ${ACC_BMAD_LINK_LIBS}")
MESSAGE("")

#-----------------------------------
# Output path definitions
#-----------------------------------
SET (LIBRARY_OUTPUT_PATH ${OUTPUT_BASEDIR}/lib)
SET (EXECUTABLE_OUTPUT_PATH ${OUTPUT_BASEDIR}/bin)
SET (CMAKE_Fortran_MODULE_DIRECTORY ${OUTPUT_BASEDIR}/modules)
SET (INCLUDE_OUTPUT_PATH ${OUTPUT_BASEDIR}/include)

#-------------------------
#   Include directories
#-------------------------
# TODO: Double each include directory entry to search for a local (../<xxx>) path and THEN
#       to a release-based path?
#
# This leaves the possibility that someone may delete the local library, and initiate a build
# that requires that library while leaving the local souce tree and include files intact.
# This new build will then perform the divergent action of linking against the release library
# but extracting constants and other header information from the LOCAL source tree.  
#
# Moved X11_INCLUDE_DIR to end of include PATH flags, as per RT#65425 and GitHub Issue #315
# for Bmad - https://github.com/DavidSagan/Bmad/issues/315

SET (MASTER_INC_DIRS
  ${INC_DIRS}
  ${OUTPUT_BASEDIR}/include
  ${OUTPUT_BASEDIR}/modules
  ${RELEASE_OUTPUT_BASEDIR}/modules
  ${RELEASE_OUTPUT_BASEDIR}/include
  ${ACC_INC_DIRS}
  ${X11_INCLUDE_DIR}
)

# If we use system HDF5 libraries, search for include directories
find_package(HDF5 COMPONENTS Fortran)
foreach(h5dir ${HDF5_Fortran_INCLUDE_DIRS})
  list(FIND CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "${h5dir}" h5found)
  if (h5found EQUAL -1)
    list(APPEND MASTER_INC_DIRS "${h5dir}")
  endif()
endforeach()

#------------------------------------------------------
# Add local include paths to search list if they exist
#------------------------------------------------------
foreach (dir ${INC_DIRS})
  STRING (FIND ${dir} "../" relative)
  STRING (REPLACE "../" "" dirnew ${dir})
  IF (${relative} EQUAL 0)
    LIST (APPEND MASTER_INC_DIRS ${dir})
    LIST (APPEND MASTER_INC_DIRS ${RELEASE_DIR}/${dirnew})
  ELSE ()
    LIST (APPEND MASTER_INC_DIRS ${dir})
  ENDIF ()
endforeach(dir)


LIST (REMOVE_DUPLICATES MASTER_INC_DIRS)
INCLUDE_DIRECTORIES (${MASTER_INC_DIRS})

#--------------------------------------------
# To avoid a Link Lib path "not found" error,
# when a Distribution Build environment is
# has been envoked and in a target project
# directory not within the BMAD_DISTRIBUTION 
# tree, create a empty lib directory.
#--------------------------------------------

IF (NOT EXISTS ${OUTPUT_BASEDIR}/lib)
  FILE (MAKE_DIRECTORY ${OUTPUT_BASEDIR}/lib)
ENDIF ()

#--------------------------------------------
# To avoid an include path "not found" error,
# when a Distribution Build environment is
# has been envoked and in a target project
# directory not within the BMAD_DISTRIBUTION 
# tree, create a empty include directory.
#--------------------------------------------

IF (NOT EXISTS ${OUTPUT_BASEDIR}/include)
  FILE (MAKE_DIRECTORY ${OUTPUT_BASEDIR}/include)
ENDIF ()

#-----------------------------------
# Link directories - order matters
# Lowest level to highest, i.e. in
# order of increasing abstraction.
#-----------------------------------

SET (MASTER_LINK_DIRS
  /usr/lib64
  ${OUTPUT_BASEDIR}/lib
  ${PACKAGES_OUTPUT_BASEDIR}/lib
  ${RELEASE_OUTPUT_BASEDIR}/lib
  ${ACC_LIB_DIRS}
)

foreach(h5dir ${HDF5_Fortran_LIBRARY_DIRS})
  list(FIND CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "${h5dir}" h5found)
  if (h5found EQUAL -1)
    list(APPEND MASTER_LINK_DIRS "${h5dir}")
  endif()
endforeach()

LINK_DIRECTORIES (${MASTER_LINK_DIRS})

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
      set_source_files_properties(${file} PROPERTIES COMPILE_FLAGS "${BASE_CXX_FLAGS} ${CPPFLAGS}")
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


set(TARGETS)


# Moved macOS Shared build specifics into a function to support SHARED only library builds, as requested in RT#63875
FUNCTION (DARWIN_SHARED)
  # This is a second refinement, may need cleanup.  Did Fix RT#47687.  -amd275 
  IF ("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
    cmake_policy (SET CMP0042 NEW)
    SET (MACOSX_RPATH TRUE)
    SET (CMAKE_SKIP_RPATH FALSE)
    SET (CMAKE_SKIP_INSTALL_RPATH FALSE)
    SET (CMAKE_MACOSX_RPATH ON)
    SET (CMAKE_INSTALL_RPATH ${RELEASE_OUTPUT_BASEDIR}/lib)
    SET (INSTALL_NAME_DIR="${RELEASE_OUTPUT_BASEDIR}/lib")
    IF (ENABLE_SHARED_ONLY) 
      SET_TARGET_PROPERTIES (${LIBNAME} PROPERTIES INSTALL_RPATH ${RELEASE_OUTPUT_BASEDIR}/lib)
      SET_TARGET_PROPERTIES (${LIBNAME} PROPERTIES MACOSX_RPATH ${RELEASE_OUTPUT_BASEDIR}/lib)
    ELSE ()
      SET_TARGET_PROPERTIES (${LIBNAME}-shared PROPERTIES INSTALL_RPATH ${RELEASE_OUTPUT_BASEDIR}/lib)
      SET_TARGET_PROPERTIES (${LIBNAME}-shared PROPERTIES MACOSX_RPATH ${RELEASE_OUTPUT_BASEDIR}/lib)
    ENDIF ()
  ENDIF ()
ENDFUNCTION ()


# Added support for SHARED only library builds, as requested in RT#63875
IF (ENABLE_SHARED_ONLY AND CREATE_SHARED)
  IF (LIBNAME)
    ADD_LIBRARY (${LIBNAME} SHARED ${sources})
    LIST (APPEND TARGETS ${LIBNAME})
    SET_TARGET_PROPERTIES (${LIBNAME} PROPERTIES OUTPUT_NAME ${LIBNAME})
    TARGET_LINK_LIBRARIES (${LIBNAME} ${SHARED_DEPS} ${SHARED_LINK_LIBS})
    DARWIN_SHARED ()
  ENDIF ()
ELSE ()
  IF (LIBNAME)
    ADD_LIBRARY (${LIBNAME} STATIC ${sources})
    LIST (APPEND TARGETS ${LIBNAME})
    SET_TARGET_PROPERTIES (${LIBNAME} PROPERTIES OUTPUT_NAME ${LIBNAME})
  ENDIF ()
ENDIF ()


#----------------------------------------------------------------
# Copy the contents of a config directory to 
# <DIR>/../config/${LIBNAME} if one exists.
#----------------------------------------------------------------
IF (IS_DIRECTORY "../config")
  message("Copying config directory contents to ${CMAKE_SOURCE_DIR}/../config/${SHORT_DIRNAME}...")
  file (MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/../config")
  EXECUTE_PROCESS (COMMAND cp -rfu ${CMAKE_SOURCE_DIR}/config/. ${CMAKE_SOURCE_DIR}/../config/${SHORT_DIRNAME})
ENDIF ()


#----------------------------------------------------------------
# For selectively producing shared object libraries (.so files).
#
# Shell variable ACC_ENABLE_SHARED needs to be set to 
#  "Y" or "true" or "1"
#    - AND - 
# set (CREATE_SHARED true) needs to be present in the individual
#  project's CMakeLists.txt file.
# 
# Now works correctly with gmake -j values greater than 1
#----------------------------------------------------------------
IF (ENABLE_SHARED)
MESSAGE ("SHARED DEPS          : ${SHARED_DEPS}")
MESSAGE ("SHARED LINKER FLAGS  : ${CMAKE_SHARED_LINKER_FLAGS}\n")
ENDIF ()
IF (ENABLE_SHARED AND CREATE_SHARED AND NOT ENABLE_SHARED_ONLY)
  ADD_LIBRARY (${LIBNAME}-shared SHARED ${sources})

  # Required when building both STATIC and SHARED Libraries are requested. 
  ADD_DEPENDENCIES (${LIBNAME}-shared ${LIBNAME}) 

  LIST (APPEND TARGETS ${LIBNAME}-shared)
  SET_TARGET_PROPERTIES (${LIBNAME}-shared PROPERTIES OUTPUT_NAME ${LIBNAME})
  TARGET_LINK_LIBRARIES (${LIBNAME}-shared ${SHARED_DEPS} ${SHARED_LINK_LIBS})
  DARWIN_SHARED ()
ENDIF ()


#---------------------------------------------------------------
# Add each TEST EXE build description file mentioned in the
# project's CMakeLists.txt file, as requested.
#---------------------------------------------------------------
IF (BUILD_TEST_EXES)
  SET (EXE_SPECS ${TEST_EXE_SPECS} ${EXE_SPECS})
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


  # Added support for SHARED only library builds, as requested in RT#63875
  IF (ENABLE_SHARED_ONLY AND CREATE_SHARED)

    foreach(dep ${DEPS})

      IF(${LIBNAME} MATCHES ${dep})
      ELSE()
          LIST(FIND TARGETS ${dep} DEP_SEEN)
          IF(${DEP_SEEN} EQUAL -1)
          IF (EXISTS ${OUTPUT_BASEDIR}/lib/lib${dep}.so)
            add_library(${dep} SHARED IMPORTED)
            LIST(APPEND TARGETS ${dep})
            set_property(TARGET ${dep} PROPERTY IMPORTED_LOCATION ${OUTPUT_BASEDIR}/lib/lib${dep}.so)
          ELSEIF (EXISTS ${PACKAGES_OUTPUT_BASEDIR}/lib/lib${dep}.so)
            add_library(${dep} SHARED IMPORTED)
            LIST(APPEND TARGETS ${dep})
            set_property(TARGET ${dep} PROPERTY IMPORTED_LOCATION ${PACKAGES_OUTPUT_BASEDIR}/lib/lib${dep}.so)
          ELSEIF (EXISTS ${RELEASE_DIR}/lib/lib${dep}.so)
            add_library(${dep} SHARED IMPORTED)
            LIST(APPEND TARGETS ${dep})
            set_property(TARGET ${dep} PROPERTY IMPORTED_LOCATION ${RELEASE_DIR}/lib/lib${dep}.so)
          ENDIF ()
      ENDIF()
    ENDIF()

    endforeach(dep)

  ELSEIF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND ENABLE_SHARED_ONLY AND CREATE_SHARED)

    foreach(dep ${DEPS})

      IF(${LIBNAME} MATCHES ${dep})
      ELSE()
          LIST(FIND TARGETS ${dep} DEP_SEEN)
          IF(${DEP_SEEN} EQUAL -1)
          IF (EXISTS ${OUTPUT_BASEDIR}/lib/lib${dep}.dylib)
            add_library(${dep} SHARED IMPORTED)
            LIST(APPEND TARGETS ${dep})
            set_property(TARGET ${dep} PROPERTY IMPORTED_LOCATION ${OUTPUT_BASEDIR}/lib/lib${dep}.dylib)
          ELSEIF (EXISTS ${PACKAGES_OUTPUT_BASEDIR}/lib/lib${dep}.dylib)
            add_library(${dep} SHARED IMPORTED)
            LIST(APPEND TARGETS ${dep})
            set_property(TARGET ${dep} PROPERTY IMPORTED_LOCATION ${PACKAGES_OUTPUT_BASEDIR}/lib/lib${dep}.dylib)
          ELSEIF (EXISTS ${RELEASE_DIR}/lib/lib${dep}.dylib)
            add_library(${dep} SHARED IMPORTED)
            LIST(APPEND TARGETS ${dep})
            set_property(TARGET ${dep} PROPERTY IMPORTED_LOCATION ${RELEASE_DIR}/lib/lib${dep}.dylib)
          ENDIF ()
  ENDIF()
      ENDIF()

    endforeach(dep)

    ELSE()

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

  ENDIF()

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
        set_source_files_properties(${file} PROPERTIES COMPILE_FLAGS "${BASE_CXX_FLAGS} ${CPPFLAGS}")
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
        set_source_files_properties(${srcfile} PROPERTIES COMPILE_FLAGS "${BASE_CXX_FLAGS} ${CFLAGS}")
      ENDIF()
      set(ext_match)
      STRING(FIND ${srcfile} ".cc" ext_match)
      IF(NOT ${ext_match} EQUAL -1)
        set_source_files_properties(${srcfile} PROPERTIES COMPILE_FLAGS "${BASE_CXX_FLAGS} ${CFLAGS}")
      ENDIF()
      set(ext_match)
      STRING(FIND ${srcfile} ".cxx" ext_match)
      IF(NOT ${ext_match} EQUAL -1)
        set_source_files_properties(${srcfile} PROPERTIES COMPILE_FLAGS "${BASE_CXX_FLAGS} ${CFLAGS}")
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

# This is a temporary attempt, needs work -amd275 
  IF ("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
    SET (CMAKE_MACOSX_RPATH ON)
    SET (CMAKE_SKIP_INSTALL_RPATH)
    SET (CMAKE_INSTALL_RPATH ${RELEASE_OUTPUT_BASEDIR}/lib)
    SET_TARGET_PROPERTIES (${EXENAME}-exe PROPERTIES INSTALL_RPATH ${RELEASE_OUTPUT_BASEDIR}/lib)
  ENDIF ()

  SET_TARGET_PROPERTIES(${EXENAME}-exe
          PROPERTIES
          OUTPUT_NAME
          ${EXENAME}
  )

  IF (DEFINED LINKER_LANGUAGE_PROP)
    SET_TARGET_PROPERTIES (${EXENAME}-exe PROPERTIES LINKER_LANGUAGE ${LINKER_LANGUAGE_PROP})
    UNSET (LINKER_LANGUAGE_PROP)
  ENDIF ()

  #----------------------------------------------
  # When linking an executable using differnet 
  # linker then the one used to build libraries 
  # in the LINK_LIBS statement, the user must 
  # set IMPLICIT_LINK_LIBS to the <lang> of the  
  # libraries in LINK_LIBS. RT#37678 
  #----------------------------------------------
  IF (DEFINED IMPLICIT_LINK_LIBS)
    SET (IMPLICIT_LINKER_LIBRARIES ${CMAKE_${IMPLICIT_LINK_LIBS}_IMPLICIT_LINK_LIBRARIES})
    MESSAGE ("Implicit ${FORTRAN_COMPILER} Linker Flags   : ${IMPLICIT_LINKER_LIBRARIES}\n")
    UNSET (IMPLICIT_LINK_LIBS)
  ENDIF ()

  # Create map file output directory if it doesn't yet exist.
  IF (IS_DIRECTORY ${OUTPUT_BASEDIR}/map)
  ELSE ()
    FILE (MAKE_DIRECTORY ${OUTPUT_BASEDIR}/map)
  ENDIF ()

  # Set up linking for the executable.
  # Always produce a map file.  It is placed in the ../<build_type>/map directory
  # created during build setup.
  set (MAPLINE "-Wl,-Map=\"${OUTPUT_BASEDIR}/map/${EXENAME}.map\"")
  IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    SET (MAPLINE "-Wl,-map -Wl,\"${OUTPUT_BASEDIR}/map/${EXENAME}.map\"")
  ENDIF ()

  SET (STATIC_FLAG "")
  SET (SHARED_FLAG "")
  IF (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    IF (FORTRAN_COMPILER MATCHES "ifort")
      # Added support for SHARED only library builds, as requested in RT#63875
      IF (ENABLE_SHARED_ONLY)
  SET (SHARED_FLAG "-Wl,-Bdynamic -fPIC")
      ELSE ()
  SET (SHARED_FLAG "-Wl,-Bdynamic -fPIC")
  SET (STATIC_FLAG "-Wl,-Bstatic -fPIC")
      ENDIF ()
    ENDIF ()

    IF ("${ACC_ENABLE_OPENMP}" AND FORTRAN_COMPILER MATCHES "gfortran")
      SET (SHARED_FLAG "${SHARED_FLAG} -fopenmp")
      SET (SHARED_DEPS "${SHARED_DEPS} gomp")
    ELSEIF () 
      SET (SHARED_FLAG "${SHARED_FLAG} -liomp5")
    ENDIF ()
  ELSEIF (${MSYS})
      # Link all libraries statically:
      SET (STATIC_FLAG "-static")
      SET (SHARED_FLAG "")
  ENDIF ()

  # Temp fix for white space issue with OPENMP Flags.
  STRING (STRIP "${SHARED_FLAG}" SHARED_FLAG)

  IF (ENABLE_SHARED AND CREATE_SHARED)
  ELSE()
    set (EXTRA_SHARED_LINK_LIBS "")
  ENDIF ()

  TARGET_LINK_LIBRARIES(${EXENAME}-exe
    ${STATIC_FLAG} ${LINK_LIBS} 
    ${SHARED_FLAG} ${SHARED_LINK_LIBS} ${EXTRA_SHARED_LINK_LIBS}
    ${X11_LIBRARIES} ${ACC_LINK_FLAGS} ${OPENMP_LINK_LIBS}
    ${LINK_FLAGS} ${MAPLINE} ${IMPLICIT_LINKER_LIBRARIES}
    )

  SET(CFLAGS)
  SET(FFLAGS)
  SET(COMPILER_FLAGS)
  SET(LINK_FLAGS)

  # Copy all header files from the project's specified include directory into the output include directory.
  SET (INCLUDE_INSTALL_DIR ${INCLUDE_OUTPUT_PATH})
  FOREACH (inc_dir ${DIR_OF_INCLUDES_TO_MOVE})
    FILE (GLOB inc_files "${inc_dir}/*.h")
    FILE (COPY ${inc_files} DESTINATION ${INCLUDE_INSTALL_DIR}) 
  ENDFOREACH (inc_dir)

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
