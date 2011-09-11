
file (GLOB IMAGEFILES "images/*.png")
file (GLOB PROFILEFILES "profiles/*.pp3")
file (GLOB LANGUAGEFILES "languages/*")
file (GLOB SOUNDFILES "sounds/*")
# THEMEDIR includes subfolders for image resources for some themes; doing the normal glob won't work.
set (THEMEDIR "themes")

if (WIN32)
    set(OPTIONSFILE "options/options.win")
elseif (APPLE)
    set(OPTIONSFILE "options/options.osx")
else (WIN32)
    set(OPTIONSFILE "options/options.lin")
endif (WIN32)

if (UNIX)
    configure_file (${CMAKE_CURRENT_SOURCE_DIR}/icons/rawtherapee.desktop.in ${CMAKE_CURRENT_BINARY_DIR}/rawtherapee.desktop)
    install (FILES ${CMAKE_CURRENT_BINARY_DIR}/rawtherapee.desktop DESTINATION ${DESKTOPDIR})
    install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/icons/hi16-app-rawtherapee.png DESTINATION ${ICONSDIR}/hicolor/16x16/apps RENAME rawtherapee.png)
    install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/icons/hi24-app-rawtherapee.png DESTINATION ${ICONSDIR}/hicolor/24x24/apps RENAME rawtherapee.png)
    install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/icons/hi32-app-rawtherapee.png DESTINATION ${ICONSDIR}/hicolor/32x32/apps RENAME rawtherapee.png)
    install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/icons/hi48-app-rawtherapee.png DESTINATION ${ICONSDIR}/hicolor/48x48/apps RENAME rawtherapee.png)
    install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/icons/hi128-app-rawtherapee.png DESTINATION ${ICONSDIR}/hicolor/128x128/apps RENAME rawtherapee.png)
    install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/icons/hi256-app-rawtherapee.png DESTINATION ${ICONSDIR}/hicolor/256x256/apps RENAME rawtherapee.png)
endif (UNIX)

install (FILES ${IMAGEFILES} DESTINATION ${DATADIR}/images)
install (FILES ${LANGUAGEFILES} DESTINATION ${DATADIR}/languages)
install (FILES ${PROFILEFILES} DESTINATION ${DATADIR}/profiles)
install (FILES ${SOUNDFILES} DESTINATION ${DATADIR}/sounds)
install (DIRECTORY ${THEMEDIR} DESTINATION ${DATADIR})
install (FILES ${OPTIONSFILE} DESTINATION ${DATADIR} PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ RENAME options)

