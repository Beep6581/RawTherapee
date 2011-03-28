# this will generate a target that will never exist, so it will (should) be executed on each build
#WARNING: Actually, only Gcc is supported

string (TOUPPER ${CMAKE_BUILD_TYPE} UPPER_CMAKE_BUILD_TYPE)

set (OUT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/AboutThisBuild.txt")

execute_process(COMMAND hg parents --template={latesttag}.{latesttagdistance} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE HG_VERSION)
execute_process(COMMAND hg parents --template={node|short} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE HG_CHANGESET)
execute_process(COMMAND hg parents --template={latesttagdistance} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE HG_TAGDISTANCE)

# build version.h from template
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/rtgui/version.h.in ${CMAKE_CURRENT_BINARY_DIR}/rtgui/version.h)

add_custom_target (AboutFile ALL
	COMMAND for /F \"tokens=*\" %%i in \('hg -R \"${CMAKE_CURRENT_SOURCE_DIR}\" branch'\) do echo Branch: %%i >${OUT_FILE}
	COMMAND echo Version: ${HG_VERSION} >>${OUT_FILE}
	COMMAND echo Changset: ${HG_CHANGESET} >>${OUT_FILE}
	COMMAND for /F \"tokens=*\" %%i in \('gcc -dumpversion'\) do echo Compiler: GCC%%i >>${OUT_FILE}
	COMMAND \(echo Processor: ${PROC_LABEL}\) >>${OUT_FILE}
	COMMAND \(echo Bit depth: ${PROC_BIT_DEPTH}\) >>${OUT_FILE}
	COMMAND \(echo Gtkmm: V${GTKMM_VERSION}\) >>${OUT_FILE}
	COMMAND \(echo Build type: ${CMAKE_BUILD_TYPE}\) >>${OUT_FILE}
	COMMAND \(echo Build flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${UPPER_CMAKE_BUILD_TYPE}} \) >>${OUT_FILE}
	COMMAND \(echo Link flags: ${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_${UPPER_CMAKE_BUILD_TYPE}} \) >>${OUT_FILE}
	COMMAND \(echo OpenMP support: ${OPTION_OMP}\) >>${OUT_FILE}
	COMMAND \(echo MMAP support: ${WITH_MYFILE_MMAP}\) >>${OUT_FILE}
	COMMAND \(echo Rawzor support: ${WITH_RAWZOR}\) >>${OUT_FILE}
	COMMENT "Creating the about file"
)
