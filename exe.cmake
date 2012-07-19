
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


add_executable(${EXENAME}-exe
    ${BUILD_EXE_TOGGLE}
    ${SRC_FILES})

SET_TARGET_PROPERTIES(${EXENAME}-exe
	PROPERTIES
	OUTPUT_NAME
	${EXENAME})

TARGET_LINK_LIBRARIES(${EXENAME}-exe
	${LINK_LIBS})