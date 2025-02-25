file(GLOB LANGUAGEFILES "languages/*")
file(GLOB SOUNDFILES "sounds/*")
file(GLOB INPUTICCFILES "iccprofiles/input/*")
file(GLOB OUTPUTICCFILES "iccprofiles/output/*")
file(GLOB DCPFILES "dcpprofiles/*")

set(PROFILESDIR "profiles")
set(THEMEDIR "themes")
set(ICONTHEMEDIR "icons")
set(IMAGESDIR "images")

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
    if(BUILD_BUNDLE)
        install(DIRECTORY "${ICONTHEMEDIR}/hicolor" DESTINATION "${ICONSDIR}")
    else()
        # For theme 'hicolor' (containing app icon), 'index.theme' shall not be installed for Linux build with 'BUILD_BUNDLE=OFF' option
        # (to avoid conflict with existing 'index.theme' distro file)
        install(DIRECTORY "${ICONTHEMEDIR}/hicolor" DESTINATION "${ICONSDIR}" PATTERN "index.theme" EXCLUDE)
    endif()
endif()

install(FILES ${LANGUAGEFILES} DESTINATION "${DATADIR}/languages")
install(FILES ${SOUNDFILES} DESTINATION "${DATADIR}/sounds")
install(FILES ${INPUTICCFILES} DESTINATION "${DATADIR}/iccprofiles/input")
install(FILES ${OUTPUTICCFILES} DESTINATION "${DATADIR}/iccprofiles/output")
install(FILES ${DCPFILES} DESTINATION "${DATADIR}/dcpprofiles")
install(FILES "cammatrices.json" DESTINATION "${DATADIR}")
install(FILES "dcraw.json" DESTINATION "${DATADIR}")
install(FILES "rt.json" DESTINATION "${DATADIR}")
install(FILES ${OPTIONSFILE} DESTINATION "${DATADIR}" PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ RENAME options)

install(DIRECTORY "${PROFILESDIR}" DESTINATION "${DATADIR}" FILES_MATCHING PATTERN "*.pp3")
install(DIRECTORY "${THEMEDIR}" DESTINATION "${DATADIR}")
install(DIRECTORY "${ICONTHEMEDIR}/rawtherapee" DESTINATION "${DATADIR}/${ICONTHEMEDIR}")
install(DIRECTORY "${IMAGESDIR}" DESTINATION "${DATADIR}")

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
