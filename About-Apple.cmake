# this will generate a target that will never exist, so it will (should) be executed on each build
#WARNING: Actually, only Gcc is supported

string (TOUPPER ${CMAKE_BUILD_TYPE} UPPER_CMAKE_BUILD_TYPE)

# we look for the hg command in this paths by order of preference
find_file(HG_CMD hg PATHS "/opt/local/bin" "/usr/local/bin" "/usr/bin")
find_file(HG_CMD hg)

if (HG_CMD STREQUAL HG_CMD-NOTFOUND)
	message(FATAL_ERROR "hg command not found!")
else (HG_CMD STREQUAL HG_CMD-NOTFOUND)
	message(STATUS "hg command found: ${HG_CMD}")
endif (HG_CMD STREQUAL HG_CMD-NOTFOUND)

set (OUT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/AboutThisBuild.txt")
set (SHELL "/bin/bash")
execute_process(COMMAND ${HG_CMD} parents --template={latesttag}.{latesttagdistance} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE HG_VERSION)
execute_process(COMMAND ${HG_CMD} parents --template={node|short} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE HG_CHANGESET)
execute_process(COMMAND ${HG_CMD} parents --template={latesttagdistance} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE HG_TAGDISTANCE)

# build version.h from template
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/rtgui/version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/rtgui/version.h)

add_custom_target (AboutFile ALL
	COMMAND rm -f ${OUT_FILE}
	COMMAND echo Branch: `${HG_CMD} -R ${CMAKE_CURRENT_SOURCE_DIR} branch` >>${OUT_FILE}
	COMMAND echo Version: ${HG_VERSION} >>${OUT_FILE}
	COMMAND echo Changeset: ${HG_CHANGESET} >>${OUT_FILE}
	COMMAND echo Compiler: GCC `gcc -dumpversion` >>${OUT_FILE}
	COMMAND echo Processor: ${PROC_LABEL} >>${OUT_FILE}
	COMMAND echo Bit depth: ${PROC_BIT_DEPTH} >>${OUT_FILE}
	COMMAND echo Gtkmm: V${GTKMM_VERSION} >>${OUT_FILE}
	COMMAND echo Build type: ${CMAKE_BUILD_TYPE} >>${OUT_FILE}
	COMMAND echo Build flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${UPPER_CMAKE_BUILD_TYPE}} >>${OUT_FILE}
	COMMAND echo Link flags: ${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_${UPPER_CMAKE_BUILD_TYPE}}  >>${OUT_FILE}
	COMMAND echo OpenMP support: ${OPTION_OMP} >>${OUT_FILE}
	COMMAND echo MMAP support: ${WITH_MYFILE_MMAP} >>${OUT_FILE}
	COMMAND echo Rawzor support: ${WITH_RAWZOR} >>${OUT_FILE}
	COMMENT "Creating the about file"
)
