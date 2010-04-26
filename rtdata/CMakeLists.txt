
file (GLOB IMAGEFILES "images/*.png")
file (GLOB PROFILEFILES "profiles/*.pp2")
file (GLOB LANGUAGEFILES "languages/*")
# THEMEDIR includes subfolders for image resources for some themes; doing the normal glob won't work.
set (THEMEDIR "themes")

if (WIN32)
    set(OPTIONSFILE "options/options.win")
elseif (APPLE)
    set(OPTIONSFILE "options/options.osx")
else (WIN32)
    set(OPTIONSFILE "options/options.lin")
endif (WIN32)


install (FILES ${IMAGEFILES} DESTINATION ${DATADIR}/images)
install (FILES ${LANGUAGEFILES} DESTINATION ${DATADIR}/languages)
install (FILES ${PROFILEFILES} DESTINATION ${DATADIR}/profiles)
install (DIRECTORY ${THEMEDIR} DESTINATION ${DATADIR})
install (FILES ${OPTIONSFILE} DESTINATION ${DATADIR} PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ RENAME options)

