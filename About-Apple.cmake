# cmakefile executed within a makefile target

# we look for the hg command in this paths by order of preference
find_file(HG_CMD hg PATHS "/opt/local/bin" "/usr/local/bin" "/usr/bin")
find_file(HG_CMD hg)

# the hg command is looked for again, at build time
if (HG_CMD STREQUAL HG_CMD-NOTFOUND)
  message(FATAL_ERROR "hg command not found!")
else (HG_CMD STREQUAL HG_CMD-NOTFOUND)
  message(STATUS "hg command found: ${HG_CMD}")
endif (HG_CMD STREQUAL HG_CMD-NOTFOUND)

set (OUT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/AboutThisBuild.txt")
set (SHELL "/bin/bash")
# there should be a better way of getting the compiler + compiler version, not restrained to Gcc
execute_process(COMMAND gcc -dumpversion OUTPUT_VARIABLE GCC_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${HG_CMD} -R ${CMAKE_CURRENT_SOURCE_DIR} branch OUTPUT_VARIABLE HG_BRANCH OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${HG_CMD} parents --template={latesttag}.{latesttagdistance} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE HG_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${HG_CMD} parents --template={node|short} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE HG_CHANGESET OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${HG_CMD} parents --template={latesttagdistance} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE HG_TAGDISTANCE OUTPUT_STRIP_TRAILING_WHITESPACE)

# build version.h from template
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/rtgui/version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/rtgui/version.h)
# build AboutThisBuild.txt from template
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/AboutThisBuild.txt.in ${CMAKE_CURRENT_SOURCE_DIR}/AboutThisBuild.txt)
