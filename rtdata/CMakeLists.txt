
file (GLOB LANGUAGEFILES "languages/*")
file (GLOB SOUNDFILES "sounds/*")
file (GLOB INPUTICCFILES "iccprofiles/input/*")
file (GLOB OUTPUTICCFILES "iccprofiles/output/*")
file (GLOB DCPFILES "dcpprofiles/*")
file (GLOB FONTS "fonts/*")
set (PROFILESDIR "profiles")
# THEMEDIR includes subfolders for image resources for some themes; doing the normal glob won't work.
set (THEMEDIR "themes")
set (IMAGESDIR "images")

if (WIN32)
    set(OPTIONSFILE "options/options.win")
elseif (APPLE)
    set(OPTIONSFILE "options/options.osx")
else (WIN32)
    set(OPTIONSFILE "options/options.lin")
endif (WIN32)

if (WIN32)
	if (CMAKE_SIZEOF_VOID_P EQUAL 4)
		set(BUILD_BIT_DEPTH 32)
		# 32 bits builds has to be installable on 64 bits system, to support WinXP/64.
		set(ARCHITECTURE_ALLOWED "x86 x64 ia64")
		# installing in 32 bits mode even on 64 bits OS and architecture
		set(INSTALL_MODE "")
		# set part of the output archive name
		set(SYSTEM_NAME "WinXP")
	elseif (CMAKE_SIZEOF_VOID_P EQUAL 8)
		set(BUILD_BIT_DEPTH 64)
		# Restricting the 64 bits builds to 64 bits systems only
		set(ARCHITECTURE_ALLOWED "x64 ia64")
		# installing in 64 bits mode for all 64 bits processors, even for itanium architecture
		set(INSTALL_MODE "x64 ia64")
		# set part of the output archive name
		set(SYSTEM_NAME "WinVista")
	endif (CMAKE_SIZEOF_VOID_P EQUAL 4)

    # If we find ReleaseInfo.cmake we use the info from there and don't need Git to be installed
    find_file(REL_INFO_FILE ReleaseInfo.cmake PATHS "${PROJECT_SOURCE_DIR}" NO_DEFAULT_PATH)
    if (REL_INFO_FILE STREQUAL REL_INFO_FILE-NOTFOUND)
        # we look for the git command in this paths by order of preference
        find_file(GIT_CMD git.exe HINTS ENV Path PATH_SUFFIXES ../)

        # Fail if Git is not installed
        if (GIT_CMD STREQUAL GIT_CMD-NOTFOUND)
            message(FATAL_ERROR "git command not found!")
        else ()
            message(STATUS "git command found: ${GIT_CMD}")
        endif ()

        execute_process(COMMAND ${GIT_CMD} -C "${PROJECT_SOURCE_DIR}" describe --tags --always OUTPUT_VARIABLE GIT_VERSION_WHOLE OUTPUT_STRIP_TRAILING_WHITESPACE)
        string(REGEX REPLACE "-.*" "" GIT_VERSION ${GIT_VERSION_WHOLE})
        string(REGEX REPLACE ".*-(.*)-g.*" "\\1" GIT_TAGDISTANCE ${GIT_VERSION_WHOLE})
    else (REL_INFO_FILE STREQUAL REL_INFO_FILE-NOTFOUND)
        include("${PROJECT_SOURCE_DIR}/ReleaseInfo.cmake")
    endif (REL_INFO_FILE STREQUAL REL_INFO_FILE-NOTFOUND)

    configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/win/InnoSetup/WindowsInnoSetup.iss.in" "${CMAKE_CURRENT_BINARY_DIR}/WindowsInnoSetup.iss")
    install (FILES "${CMAKE_CURRENT_BINARY_DIR}/WindowsInnoSetup.iss" DESTINATION ${BINDIR})
endif (WIN32)

if (UNIX)
    configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/icons/rawtherapee.desktop.in" "${CMAKE_CURRENT_BINARY_DIR}/rawtherapee.desktop")
    install (FILES "${CMAKE_CURRENT_BINARY_DIR}/rawtherapee.desktop" DESTINATION ${DESKTOPDIR})
    install (FILES "${CMAKE_CURRENT_SOURCE_DIR}/icons/hi16-app-rawtherapee.png" DESTINATION "${ICONSDIR}/hicolor/16x16/apps" RENAME rawtherapee.png)
    install (FILES "${CMAKE_CURRENT_SOURCE_DIR}/icons/hi24-app-rawtherapee.png" DESTINATION "${ICONSDIR}/hicolor/24x24/apps" RENAME rawtherapee.png)
    install (FILES "${CMAKE_CURRENT_SOURCE_DIR}/icons/hi32-app-rawtherapee.png" DESTINATION "${ICONSDIR}/hicolor/32x32/apps" RENAME rawtherapee.png)
    install (FILES "${CMAKE_CURRENT_SOURCE_DIR}/icons/hi48-app-rawtherapee.png" DESTINATION "${ICONSDIR}/hicolor/48x48/apps" RENAME rawtherapee.png)
    install (FILES "${CMAKE_CURRENT_SOURCE_DIR}/icons/hi128-app-rawtherapee.png" DESTINATION "${ICONSDIR}/hicolor/128x128/apps" RENAME rawtherapee.png)
    install (FILES "${CMAKE_CURRENT_SOURCE_DIR}/icons/hi256-app-rawtherapee.png" DESTINATION "${ICONSDIR}/hicolor/256x256/apps" RENAME rawtherapee.png)
endif (UNIX)

install (FILES ${IMAGEFILES} DESTINATION "${DATADIR}/images")
install (FILES ${LANGUAGEFILES} DESTINATION "${DATADIR}/languages")
install (FILES ${SOUNDFILES} DESTINATION "${DATADIR}/sounds")
install (FILES ${INPUTICCFILES} DESTINATION "${DATADIR}/iccprofiles/input")
install (FILES ${OUTPUTICCFILES} DESTINATION "${DATADIR}/iccprofiles/output")
install (FILES ${DCPFILES} DESTINATION "${DATADIR}/dcpprofiles")
if (WIN32)
    install (FILES ${FONTS} DESTINATION "${DATADIR}/fonts")
endif (WIN32)
install (DIRECTORY ${PROFILESDIR} DESTINATION "${DATADIR}" FILES_MATCHING PATTERN "*.pp3")
install (DIRECTORY ${THEMEDIR} DESTINATION "${DATADIR}")
install (DIRECTORY ${IMAGESDIR} DESTINATION "${DATADIR}" FILES_MATCHING PATTERN "index.theme")
install (DIRECTORY ${IMAGESDIR} DESTINATION "${DATADIR}" FILES_MATCHING PATTERN "*.png")
install (FILES ${OPTIONSFILE} DESTINATION "${DATADIR}" PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ RENAME options)

if (WIN32)
  install (FILES "${PROJECT_SOURCE_DIR}/doc/manpage/rawtherapee.1" DESTINATION "${CMAKE_INSTALL_PREFIX}/share/man/man1")
endif (WIN32)

if (APPLE)
    # CMake escapes first item quote character. Do not remove 'DUMMY_VARIABLE='.
    set (MACOSX_BUNDLE_COMMAND  DUMMY_VARIABLE=
                                PROJECT_NAME="${PROJECT_NAME}"
                                PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
                                CMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}"
                                GTK_PREFIX="${GTK_PREFIX}")
    if (PROC_BIT_DEPTH MATCHES 32)
        list (APPEND MACOSX_BUNDLE_COMMAND PROC_BIT_DEPTH=32)
    elseif (PROC_BIT_DEPTH MATCHES 64)
        list (APPEND MACOSX_BUNDLE_COMMAND PROC_BIT_DEPTH=64)
    endif (PROC_BIT_DEPTH MATCHES 32)
    list (APPEND MACOSX_BUNDLE_COMMAND sh "${PROJECT_SOURCE_DIR}/rtdata/osx/macosx_bundle.sh")
    
    add_custom_target(macosx_bundle
        COMMAND ${MACOSX_BUNDLE_COMMAND}
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
        COMMENT "Creating Mac OS X bundle")
endif (APPLE)
