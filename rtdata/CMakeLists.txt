
file (GLOB IMAGEFILES "images/*.png")
file (GLOB PROFILEFILES "profiles/*.pp2")
# unfortunatelly GLOB * matches directories too which breas cmake build from
# sources handled by subversion, so better stick with manual list for these
#file(GLOB LANGUAGEFILES "languages/*")
#file(GLOB THEMEFILES "themes/*")
set (LANGUAGEFILES 
	"languages/Chinese (Simplified)"
	"languages/Chinese (Traditional)"
	"languages/English (UK)"
	"languages/English (US)"
	"languages/Norsk BM"
	"languages/Polish (Latin Characters)"
	"languages/Portugues (Brasil)")
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

