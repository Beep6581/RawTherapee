
file (GLOB IMAGEFILES "images/*.png")
file (GLOB PROFILEFILES "profiles/*.pp2")
# unfortunatelly GLOB * matches directories too which breas cmake build from
# sources handled by subversion, so better stick with manual list for these
#file(GLOB LANGUAGEFILES "languages/*")
#file(GLOB THEMEFILES "themes/*")
set (LANGUAGEFILES 
	"languages/Catala"
	"languages/Chinese (Simplified)"
	"languages/Chinese (Traditional)"
	"languages/Czech"
	"languages/Dansk"
	"languages/Deutsch"
	"languages/English (UK)"
	"languages/English (US)"
	"languages/Espanol"
	"languages/Euskara"
	"languages/Francais"
	"languages/Greek"
	"languages/Hebrew"
	"languages/Italian"
	"languages/Japanese"
	"languages/Latvian"
	"languages/Magyar"
	"languages/Nederlands"
	"languages/Norsk BM"
	"languages/Polish"
	"languages/Polish (Latin Characters)"
	"languages/Portugues (Brasil)"
	"languages/Russian"
	"languages/Slovak"
	"languages/Suomi"
	"languages/Swedish"
	"languages/Turkish")
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

