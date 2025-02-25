# cmakefile executed within a makefile target

if(APPLE)
    set(PROJECT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/..")
endif()

# If we find ReleaseInfo.cmake we use the info from there and don't need Git to be installed
find_file(REL_INFO_FILE ReleaseInfo.cmake PATHS "${PROJECT_SOURCE_DIR}" NO_DEFAULT_PATH)
if(REL_INFO_FILE STREQUAL REL_INFO_FILE-NOTFOUND)
    # we look for the git command in this paths by order of preference
    if(WIN32)
        find_program(GIT_CMD git.exe HINTS ENV Path PATH_SUFFIXES ../)
    elseif(APPLE)
        find_program(GIT_CMD git PATHS "/opt/local/bin" "/usr/local/bin" "/usr/bin")
        find_program(GIT_CMD git)
        set(SHELL "/bin/bash")
    else() # Linux
        find_program(GIT_CMD git)
        set(SHELL "/bin/bash")
    endif()

    # Fail if Git is not installed
    if(GIT_CMD STREQUAL GIT_CMD-NOTFOUND)
        message(FATAL_ERROR "git command not found!")
    else()
        message(STATUS "git command found: ${GIT_CMD}")
    endif()

    # Get version description.
    # Depending on whether you checked out a branch (dev) or a tag (release),
    # "git describe" will return "5.0-gtk2-2-g12345678" or "5.0-gtk2", respectively.
    execute_process(COMMAND ${GIT_CMD} describe --tags --always OUTPUT_VARIABLE GIT_DESCRIBE OUTPUT_STRIP_TRAILING_WHITESPACE WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")

    # Get branch name.
    # Will return empty if you checked out a commit or tag. Empty string handled later.
    execute_process(COMMAND ${GIT_CMD} symbolic-ref --short -q HEAD OUTPUT_VARIABLE GIT_BRANCH OUTPUT_STRIP_TRAILING_WHITESPACE WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")

    # Get commit hash.
    execute_process(COMMAND ${GIT_CMD} rev-parse --short --verify HEAD OUTPUT_VARIABLE GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")

    # Get commit date, YYYY-MM-DD.
    execute_process(COMMAND ${GIT_CMD} show -s --format=%cd --date=format:%Y-%m-%d OUTPUT_VARIABLE GIT_COMMIT_DATE OUTPUT_STRIP_TRAILING_WHITESPACE WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")

    # Get number of commits since tagging. This is what "GIT_DESCRIBE" uses.
    # Works when checking out branch, tag or commit.
    # Get a list of all tags in repo:
    execute_process(COMMAND ${GIT_CMD} tag --merged HEAD OUTPUT_VARIABLE GIT_TAG WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
    # Replace newlines with semicolons so that it can be split:
    string(REPLACE "\n" ";" GIT_TAG_LIST "${GIT_TAG}")
    execute_process(COMMAND ${GIT_CMD} rev-list --count HEAD --not ${GIT_TAG_LIST} OUTPUT_VARIABLE GIT_COMMITS_SINCE_TAG OUTPUT_STRIP_TRAILING_WHITESPACE WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")

    # Get number of commits since branching.
    # Works when checking out branch, tag or commit.
    execute_process(COMMAND ${GIT_CMD} rev-list --count HEAD --not --tags OUTPUT_VARIABLE GIT_COMMITS_SINCE_BRANCH OUTPUT_STRIP_TRAILING_WHITESPACE WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")

    # If user checked-out something which is not a branch, use the description as branch name.
    if(GIT_BRANCH STREQUAL "")
        set(GIT_BRANCH "${GIT_DESCRIBE}")
    endif()

    # Create numeric version.
    # This version is nonsense, either don't use it at all or use it only where you have no other choice, e.g. Inno Setup's VersionInfoVersion.
    # Strip everything after hyphen, e.g. "5.0-gtk2" -> "5.0", "5.1-rc1" -> "5.1" (ergo BS).
    if(GIT_COMMITS_SINCE_TAG STREQUAL "")
        set(GIT_NUMERIC_VERSION_BS "0.0.0")
    else()
        string(REGEX REPLACE "-.*" "" GIT_NUMERIC_VERSION_BS ${GIT_DESCRIBE})
        set(GIT_NUMERIC_VERSION_BS "${GIT_NUMERIC_VERSION_BS}.${GIT_COMMITS_SINCE_TAG}")
    endif()

    execute_process(COMMAND uname -mrs OUTPUT_VARIABLE BUILDINFO_OS OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND date -Ru OUTPUT_VARIABLE BUILDINFO_DATE OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND date +%s OUTPUT_VARIABLE BUILDINFO_EPOCH OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND uuidgen COMMAND tr "A-Z" "a-z" OUTPUT_VARIABLE BUILDINFO_UUID OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "Git checkout information:")
    message(STATUS "    Commit description:    ${GIT_DESCRIBE}")
    message(STATUS "    Branch:                ${GIT_BRANCH}")
    message(STATUS "    Commit:                ${GIT_COMMIT}")
    message(STATUS "    Commit date:           ${GIT_COMMIT_DATE}")
    message(STATUS "    Commits since tag:     ${GIT_COMMITS_SINCE_TAG}")
    message(STATUS "    Commits since branch:  ${GIT_COMMITS_SINCE_BRANCH}")
    message(STATUS "    Version (unreliable):  ${GIT_NUMERIC_VERSION_BS}")
    message(STATUS "Build information:")
    message(STATUS "    Build OS:              ${BUILDINFO_OS}")
    message(STATUS "    Build date:            ${BUILDINFO_DATE} UTC")
    message(STATUS "    Epoch:                 ${BUILDINFO_EPOCH}")
    message(STATUS "    UUID:                  ${BUILDINFO_UUID}")

    if(NOT DEFINED CACHE_NAME_SUFFIX)
        set(CACHE_NAME_SUFFIX "${GIT_DESCRIBE}")
        message(STATUS "CACHE_NAME_SUFFIX was not defined, it is now \"${CACHE_NAME_SUFFIX}\"")
    else()
        message(STATUS "CACHE_NAME_SUFFIX is \"${CACHE_NAME_SUFFIX}\"")
    endif()
else()
    include("${PROJECT_SOURCE_DIR}/ReleaseInfo.cmake")
endif()

if(WIN32)
    if(BIT_DEPTH EQUAL 4)
        set(BUILD_BIT_DEPTH 32)
        # 32 bits builds has to be installable on 64 bits system, to support WinXP/64.
        set(ARCHITECTURE_ALLOWED "x86 x64")
        # installing in 32 bits mode even on 64 bits OS and architecture
        set(INSTALL_MODE "")
    elseif(BIT_DEPTH EQUAL 8)
        set(BUILD_BIT_DEPTH 64)
        # Restricting the 64 bits builds to 64 bits systems only
        set(ARCHITECTURE_ALLOWED "x64 arm64")
        # installing in 64 bits mode for all 64 bits processors, even for itanium architecture
        set(INSTALL_MODE "x64 arm64")
    endif()
    # set part of the output archive name
    set(SYSTEM_NAME "WinVista")

    configure_file("${PROJECT_SOURCE_DIR}/tools/win/InnoSetup/WindowsInnoSetup.iss.in" "${CMAKE_BINARY_DIR}/rtdata/WindowsInnoSetup.iss")
endif()

# build version.h from template
configure_file("${PROJECT_SOURCE_DIR}/rtgui/version.h.in" "${CMAKE_BINARY_DIR}/rtgui/version.h")
# build AboutThisBuild.txt from template
configure_file("${PROJECT_SOURCE_DIR}/AboutThisBuild.txt.in" "${CMAKE_BINARY_DIR}/AboutThisBuild.txt")
