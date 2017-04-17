
file(GLOB LANGUAGEFILES "languages/*")
file(GLOB SOUNDFILES "sounds/*")
file(GLOB INPUTICCFILES "iccprofiles/input/*")
file(GLOB OUTPUTICCFILES "iccprofiles/output/*")
file(GLOB DCPFILES "dcpprofiles/*")
file(GLOB FONTS "fonts/*")

set(PROFILESDIR "profiles")
set(IMAGESDIR "images")

# THEMEDIR includes subfolders for image resources for some themes; doing the normal glob won't work.
set(THEMEDIR "themes")

if(WIN32)
    set(OPTIONSFILE "options/options.win")
elseif(APPLE)
    set(OPTIONSFILE "options/options.osx")
else()
    set(OPTIONSFILE "options/options.lin")
endif()

if(WIN32)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/WindowsInnoSetup.iss" DESTINATION ${BINDIR})
endif()

if(UNIX)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/icons/rawtherapee.desktop.in" "${CMAKE_CURRENT_BINARY_DIR}/rawtherapee.desktop")
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/rawtherapee.desktop" DESTINATION ${DESKTOPDIR})
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/icons/hi16-app-rawtherapee.png" DESTINATION "${ICONSDIR}/hicolor/16x16/apps" RENAME rawtherapee.png)
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/icons/hi24-app-rawtherapee.png" DESTINATION "${ICONSDIR}/hicolor/24x24/apps" RENAME rawtherapee.png)
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/icons/hi32-app-rawtherapee.png" DESTINATION "${ICONSDIR}/hicolor/32x32/apps" RENAME rawtherapee.png)
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/icons/hi48-app-rawtherapee.png" DESTINATION "${ICONSDIR}/hicolor/48x48/apps" RENAME rawtherapee.png)
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/icons/hi128-app-rawtherapee.png" DESTINATION "${ICONSDIR}/hicolor/128x128/apps" RENAME rawtherapee.png)
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/icons/hi256-app-rawtherapee.png" DESTINATION "${ICONSDIR}/hicolor/256x256/apps" RENAME rawtherapee.png)
endif()

install(FILES ${IMAGEFILES} DESTINATION "${DATADIR}/images")
install(FILES ${LANGUAGEFILES} DESTINATION "${DATADIR}/languages")
install(FILES ${SOUNDFILES} DESTINATION "${DATADIR}/sounds")
install(FILES ${INPUTICCFILES} DESTINATION "${DATADIR}/iccprofiles/input")
install(FILES ${OUTPUTICCFILES} DESTINATION "${DATADIR}/iccprofiles/output")
install(FILES ${DCPFILES} DESTINATION "${DATADIR}/dcpprofiles")
install(FILES ${OPTIONSFILE} DESTINATION "${DATADIR}" PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ RENAME options)
if(WIN32)
    install(FILES ${FONTS} DESTINATION "${DATADIR}/fonts")
endif()

install(DIRECTORY ${PROFILESDIR} DESTINATION "${DATADIR}" FILES_MATCHING PATTERN "*.pp3")
install(DIRECTORY ${THEMEDIR} DESTINATION "${DATADIR}")
install(DIRECTORY ${IMAGESDIR} DESTINATION "${DATADIR}" FILES_MATCHING PATTERN "index.theme")
install(DIRECTORY ${IMAGESDIR} DESTINATION "${DATADIR}" FILES_MATCHING PATTERN "*.png")

if(APPLE)
    # CMake escapes first item quote character. Do not remove 'DUMMY_VARIABLE='
    set(MACOSX_BUNDLE_COMMAND DUMMY_VARIABLE=
        PROJECT_NAME="${PROJECT_NAME}"
        PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
        CMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}"
        GTK_PREFIX="${GTK_PREFIX}")
    if(PROC_BIT_DEPTH MATCHES 32)
        list(APPEND MACOSX_BUNDLE_COMMAND PROC_BIT_DEPTH=32)
    else()
        list(APPEND MACOSX_BUNDLE_COMMAND PROC_BIT_DEPTH=64)
    endif()
    list(APPEND MACOSX_BUNDLE_COMMAND sh "${PROJECT_SOURCE_DIR}/tools/osx/macosx_bundle.sh")

    add_custom_target(macosx_bundle
        COMMAND ${MACOSX_BUNDLE_COMMAND}
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
        COMMENT "Creating macOS bundle")
endif()
