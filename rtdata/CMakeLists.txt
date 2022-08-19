file(GLOB LANGUAGEFILES "languages/*")
file(GLOB SOUNDFILES "sounds/*")
file(GLOB INPUTICCFILES "iccprofiles/input/*")
file(GLOB OUTPUTICCFILES "iccprofiles/output/*")
file(GLOB DCPFILES "dcpprofiles/*")
file(GLOB FONTS "fonts/*")

set(PROFILESDIR "profiles")
set(THEMEDIR "themes")

file(GLOB IMG_SVG LIST_DIRECTORIES false "images/svg/*.svg")
file(GLOB IMG_ICO LIST_DIRECTORIES false "images/*")

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
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/rawtherapee.desktop.in" "${CMAKE_CURRENT_BINARY_DIR}/rawtherapee.desktop")
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/rawtherapee.desktop" DESTINATION ${DESKTOPDIR})
endif()

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/images/png/rawtherapee-logo-16.png" DESTINATION "${ICONSDIR}/rawtherapee/16x16/app_icons" RENAME rawtherapee.png)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/images/png/rawtherapee-logo-24.png" DESTINATION "${ICONSDIR}/rawtherapee/24x24/app_icons" RENAME rawtherapee.png)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/images/png/rawtherapee-logo-32.png" DESTINATION "${ICONSDIR}/rawtherapee/32x32/app_icons" RENAME rawtherapee.png)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/images/png/rawtherapee-logo-48.png" DESTINATION "${ICONSDIR}/rawtherapee/48x48/app_icons" RENAME rawtherapee.png)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/images/png/rawtherapee-logo-64.png" DESTINATION "${ICONSDIR}/rawtherapee/64x64/app_icons" RENAME rawtherapee.png)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/images/png/rawtherapee-logo-128.png" DESTINATION "${ICONSDIR}/rawtherapee/128x128/app_icons" RENAME rawtherapee.png)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/images/png/rawtherapee-logo-256.png" DESTINATION "${ICONSDIR}/rawtherapee/256x256/app_icons" RENAME rawtherapee.png)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/images/rt-logo.svg" DESTINATION "${ICONSDIR}/rawtherapee/scalable/app_icons" RENAME rawtherapee.svg)

install(FILES ${LANGUAGEFILES} DESTINATION "${DATADIR}/languages")
install(FILES ${SOUNDFILES} DESTINATION "${DATADIR}/sounds")
install(FILES ${INPUTICCFILES} DESTINATION "${DATADIR}/iccprofiles/input")
install(FILES ${OUTPUTICCFILES} DESTINATION "${DATADIR}/iccprofiles/output")
install(FILES ${DCPFILES} DESTINATION "${DATADIR}/dcpprofiles")
install(FILES ${OPTIONSFILE} DESTINATION "${DATADIR}" PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ RENAME options)
if(WIN32)
    install(FILES ${FONTS} DESTINATION "${DATADIR}/fonts")
endif()

install(DIRECTORY "${PROFILESDIR}" DESTINATION "${DATADIR}" FILES_MATCHING PATTERN "*.pp3")
install(DIRECTORY "${THEMEDIR}" DESTINATION "${DATADIR}")

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/images/png/empty.png" DESTINATION "${ICONSDIR}/rawtherapee/24x24/apps")
install(FILES ${IMG_SVG} DESTINATION "${ICONSDIR}/rawtherapee/scalable/apps")
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/images/svg/Index.theme" DESTINATION "${ICONSDIR}/rawtherapee")
install(FILES ${IMG_ICO} DESTINATION "${DATADIR}/images")

if(APPLE)
    # CMake escapes first item quote character. Do not remove 'DUMMY_VARIABLE='
    set(MACOSX_BUNDLE_COMMAND DUMMY_VARIABLE=
        PROJECT_NAME=${PROJECT_NAME}
        PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
        CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        GTK_PREFIX="${GTK_PREFIX}")
    list(APPEND MACOSX_BUNDLE_COMMAND PROC_BIT_DEPTH=64)
    list(APPEND MACOSX_BUNDLE_COMMAND sh ${PROJECT_SOURCE_DIR}/tools/osx/macosx_bundle.sh)

    add_custom_target(macosx_bundle
        COMMAND ${MACOSX_BUNDLE_COMMAND}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/build
        COMMENT "Creating macOS bundle")
endif()
