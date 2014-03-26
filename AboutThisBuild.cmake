# cmakefile executed within a makefile target

# If we find ReleaseInfo.cmake we use the info from there and don't need Mercurial to be installed
find_file(REL_INFO_FILE ReleaseInfo.cmake PATHS "${PROJECT_SOURCE_DIR}" NO_DEFAULT_PATH)
if (REL_INFO_FILE STREQUAL REL_INFO_FILE-NOTFOUND)
  # we look for the hg command in this paths by order of preference
  if (WIN32)
    find_file(HG_CMD hg.exe HINTS ENV Path PATH_SUFFIXES ../)
  elseif (APPLE)
    find_file(HG_CMD hg PATHS "/opt/local/bin" "/usr/local/bin" "/usr/bin")
    find_file(HG_CMD hg)
    set (SHELL "/bin/bash")
  else (WIN32) # Linux
    find_file(HG_CMD hg)
    set (SHELL "/bin/bash")
  endif (WIN32)

  # Fail if Mercurial is not installed
  if (HG_CMD STREQUAL HG_CMD-NOTFOUND)
    message(FATAL_ERROR "hg command not found!")
  else (HG_CMD STREQUAL HG_CMD-NOTFOUND)
    message(STATUS "hg command found: ${HG_CMD}")
  endif (HG_CMD STREQUAL HG_CMD-NOTFOUND)

  execute_process(COMMAND ${HG_CMD} -R "${PROJECT_SOURCE_DIR}" branch OUTPUT_VARIABLE HG_BRANCH OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${HG_CMD} -R "${PROJECT_SOURCE_DIR}" parents --template={latesttag}.{latesttagdistance} WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" OUTPUT_VARIABLE HG_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${HG_CMD} -R "${PROJECT_SOURCE_DIR}" parents --template={node|short} WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" OUTPUT_VARIABLE HG_CHANGESET OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${HG_CMD} -R "${PROJECT_SOURCE_DIR}" parents --template={latesttagdistance} WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" OUTPUT_VARIABLE HG_TAGDISTANCE OUTPUT_STRIP_TRAILING_WHITESPACE)
  if (CACHE_NAME_SUFFIX STREQUAL "")
    execute_process(COMMAND ${HG_CMD} -R "${PROJECT_SOURCE_DIR}" parents --template={latesttag} WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" OUTPUT_VARIABLE CACHE_NAME_SUFFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif (CACHE_NAME_SUFFIX STREQUAL "")
else (REL_INFO_FILE STREQUAL REL_INFO_FILE-NOTFOUND)
  include("${PROJECT_SOURCE_DIR}/ReleaseInfo.cmake")
endif (REL_INFO_FILE STREQUAL REL_INFO_FILE-NOTFOUND)

if (VERSION_SUFFIX STREQUAL "")
  set (HG_VERSION_SUFFIX "${HG_VERSION}")
else (VERSION_SUFFIX STREQUAL "")
  set (HG_VERSION_SUFFIX "${HG_VERSION} ${VERSION_SUFFIX}")
endif ()

# build version.h from template
configure_file ("${PROJECT_SOURCE_DIR}/rtgui/version.h.in" "${CMAKE_BINARY_DIR}/rtgui/version.h")
# build AboutThisBuild.txt from template
configure_file ("${PROJECT_SOURCE_DIR}/AboutThisBuild.txt.in" "${CMAKE_BINARY_DIR}/AboutThisBuild.txt")
