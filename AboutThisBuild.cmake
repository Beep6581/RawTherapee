# cmakefile executed within a makefile target

# If we find ReleaseInfo.cmake we use the info from there and don't need Git to be installed
find_file(REL_INFO_FILE ReleaseInfo.cmake PATHS "${PROJECT_SOURCE_DIR}" NO_DEFAULT_PATH)
if (REL_INFO_FILE STREQUAL REL_INFO_FILE-NOTFOUND)
    # we look for the git command in this paths by order of preference
    if (WIN32)
        find_file(GIT_CMD git.exe HINTS ENV Path PATH_SUFFIXES ../)
    elseif (APPLE)
        find_file(GIT_CMD git PATHS "/opt/local/bin" "/usr/local/bin" "/usr/bin")
        find_file(GIT_CMD git)
        set (SHELL "/bin/bash")
    else (WIN32) # Linux
        find_file(GIT_CMD git)
        set (SHELL "/bin/bash")
    endif (WIN32)

    # Fail if Git is not installed
    if (GIT_CMD STREQUAL GIT_CMD-NOTFOUND)
        message(FATAL_ERROR "git command not found!")
    else ()
        message(STATUS "git command found: ${GIT_CMD}")
    endif ()

    execute_process(COMMAND ${GIT_CMD} symbolic-ref --short -q HEAD OUTPUT_VARIABLE GIT_BRANCH OUTPUT_STRIP_TRAILING_WHITESPACE WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
    execute_process(COMMAND ${GIT_CMD} describe --tags --always OUTPUT_VARIABLE GIT_VERSION_WHOLE OUTPUT_STRIP_TRAILING_WHITESPACE WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
    string(REGEX REPLACE "-g.*" "" GIT_VERSION ${GIT_VERSION_WHOLE})
    string(REPLACE "-" "." GIT_VERSION ${GIT_VERSION})
    execute_process(COMMAND ${GIT_CMD} rev-parse --verify HEAD OUTPUT_VARIABLE GIT_CHANGESET OUTPUT_STRIP_TRAILING_WHITESPACE WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
    string(REGEX REPLACE ".*-(.*)-g.*" "\\1" GIT_TAGDISTANCE ${GIT_VERSION_WHOLE})
    if (NOT DEFINED CACHE_NAME_SUFFIX)
        string(REGEX REPLACE "-.*" "" CACHE_NAME_SUFFIX ${GIT_VERSION_WHOLE})
        message(STATUS "CACHE_NAME_SUFFIX was not defined, it is now \"${CACHE_NAME_SUFFIX}\"")
    elseif (CACHE_NAME_SUFFIX STREQUAL "latesttag")
        string(REGEX REPLACE "-.*" "" CACHE_NAME_SUFFIX ${GIT_VERSION_WHOLE})
        message(STATUS "CACHE_NAME_SUFFIX was \"latesttag\", it is now \"${CACHE_NAME_SUFFIX}\"")
    else ()
        message(STATUS "CACHE_NAME_SUFFIX is \"${CACHE_NAME_SUFFIX}\"")
    endif ()
else (REL_INFO_FILE STREQUAL REL_INFO_FILE-NOTFOUND)
    include("${PROJECT_SOURCE_DIR}/ReleaseInfo.cmake")
endif (REL_INFO_FILE STREQUAL REL_INFO_FILE-NOTFOUND)

if (VERSION_SUFFIX STREQUAL "")
    set (GIT_VERSION_SUFFIX "${GIT_VERSION}")
else ()
    set (GIT_VERSION_SUFFIX "${GIT_VERSION} ${VERSION_SUFFIX}")
endif ()

# build version.h from template
configure_file ("${PROJECT_SOURCE_DIR}/rtgui/version.h.in" "${CMAKE_BINARY_DIR}/rtgui/version.h")
# build AboutThisBuild.txt from template
configure_file ("${PROJECT_SOURCE_DIR}/AboutThisBuild.txt.in" "${CMAKE_BINARY_DIR}/AboutThisBuild.txt")
