
find_package(PkgConfig)

pkg_check_modules (GLIB2   REQUIRED glib-2.0>=2.16)
pkg_check_modules (GLIBMM  REQUIRED glibmm-2.4>=2.16)
pkg_check_modules (GTK     REQUIRED gtk+-2.0>=2.12)
pkg_check_modules (GTKMM   REQUIRED gtkmm-2.4>=2.12)
pkg_check_modules (GIO     REQUIRED gio-2.0>=2.16)
pkg_check_modules (GIOMM   REQUIRED giomm-2.4>=2.12)

SET (BASESOURCEFILES
    batchtoolpanelcoord.cc paramsedited.cc cropwindow.cc previewhandler.cc previewwindow.cc navigator.cc indclippedpanel.cc filterpanel.cc
    cursormanager.cc rtwindow.cc renamedlg.cc recentbrowser.cc placesbrowser.cc filepanel.cc editorpanel.cc batchqueuepanel.cc
    ilabel.cc thumbbrowserbase.cc adjuster.cc filebrowserentry.cc filebrowser.cc filethumbnailbuttonset.cc
    cachemanager.cc cacheimagedata.cc
    clipboard.cc thumbimageupdater.cc bqentryupdater.cc
    coarsepanel.cc cacorrection.cc colorshift.cc hlrec.cc chmixer.cc
    colorboost.cc resize.cc icmpanel.cc crop.cc shadowshighlights.cc
    colordenoise.cc
    exifpanel.cc
    sharpening.cc
    whitebalance.cc   vignetting.cc rotate.cc distortion.cc
    crophandler.cc curveeditor.cc dirbrowser.cc
    filecatalog.cc
    histogrampanel.cc history.cc  imagearea.cc
    imageareapanel.cc iptcpanel.cc lcurve.cc lumadenoise.cc main.cc
    multilangmgr.cc mycurve.cc options.cc
    preferences.cc profilepanel.cc saveasdlg.cc
    saveformatpanel.cc splash.cc
    thumbnail.cc tonecurve.cc toolbar.cc
    guiutils.cc zoompanel.cc toolpanelcoord.cc
    thumbbrowserentrybase.cc batchqueueentry.cc
    batchqueue.cc lwbutton.cc lwbuttonset.cc
    batchqueuebuttonset.cc browserfilter.cc exiffiltersettings.cc
    profilestore.cc partialpastedlg.cc)

IF (WIN32)
    SET (EXTRA_LIBDIR "${CMAKE_CURRENT_SOURCE_DIR}/../lib")
    SET (EXTRA_INCDIR "${CMAKE_CURRENT_SOURCE_DIR}/../winclude")
    SET (EXTRA_SRC "windirmonitor.cc myicon.o")
    SET (EXTRA_LIB "ws2_32")

include_directories (/usr/local/lib ../rtengine . ../rtexif ${EXTRA_INCDIR} ${GLIB2_INCLUDE_DIRS} ${GLIBMM_INCLUDE_DIRS}
    ${GTK_INCLUDE_DIRS} ${GTKMM_INCLUDE_DIRS} ${GIO_INCLUDE_DIRS} ${GIOMM_INCLUDE_DIRS})
link_directories (. ../rtexif ${EXTRA_LIBDIR} ${GLIB2_LIBRARY_DIRS} ${GLIBMM_LIBRARY_DIRS}
    ${GTK_LIBRARY_DIRS} ${GTKMM_LIBRARY_DIRS} ${GIO_LIBRARY_DIRS} ${GIOMM_LIBRARY_DIRS})
#set_target_properties (rth PROPERTIES LINK_FLAGS "-mwindows")

    add_executable (rth windirmonitor.cc myicon.o ${BASESOURCEFILES})

    set(IPTCDATA_LIBRARIES iptcdata)
    set(LCMS_LIBRARIES liblcms.a)
    set(GTHREAD_LIBRARIES gthread-2.0)
    set(GOBJECTLIBRARIES gobject-2.0)
    set(JPEG_LIBRARIES libjpeg.a)
    set(ZLIB_LIBRARIES libz.a)
    set(TIFF_LIBRARIES libtiff.a)
    set(PNG_LIBRARIES libpng.a)

ELSE (WIN32)
    pkg_check_modules (IPTCDATA REQUIRED libiptcdata)
    pkg_check_modules (LCMS REQUIRED lcms)
    pkg_check_modules (PNG REQUIRED libpng)
    pkg_check_modules (GTHREAD REQUIRED gthread-2.0>=2.16)
    pkg_check_modules (GOBJECT REQUIRED gobject-2.0>=2.16)
    find_library(JPEG_LIBRARIES NAMES jpeg)
    if (NOT JPEG_LIBRARIES)
        message(SEND_ERROR "libjpeg was not found!")
    endif (NOT JPEG_LIBRARIES)
    find_library(TIFF_LIBRARIES NAMES tiff)
    if (NOT TIFF_LIBRARIES)
        message(SEND_ERROR "libtiff was not found!")
    endif (NOT TIFF_LIBRARIES)
    find_library(ZLIB_LIBRARIES NAMES z)
    if (NOT ZLIB_LIBRARIES)
        message(SEND_ERROR "libz was not found!")
    endif (NOT ZLIB_LIBRARIES)
    if (WITH_RAWZOR)
    IF (CMAKE_SIZEOF_VOID_P EQUAL 4)
        SET (EXTRA_INCDIR "${CMAKE_CURRENT_SOURCE_DIR}/../rawzor_lin32")
	    SET (EXTRA_LIBDIR "${CMAKE_CURRENT_SOURCE_DIR}/../rawzor_lin32")
    ELSEIF (CMAKE_SIZEOF_VOID_P EQUAL 8)
        SET (EXTRA_INCDIR "${CMAKE_CURRENT_SOURCE_DIR}/../rawzor_lin64")
	    SET (EXTRA_LIBDIR "${CMAKE_CURRENT_SOURCE_DIR}/../rawzor_lin64")
    ENDIF (CMAKE_SIZEOF_VOID_P EQUAL 4)
    endif (WITH_RAWZOR)

include_directories (/usr/local/lib ${CMAKE_CURRENT_SOURCE_DIR}/../rtengine . ${CMAKE_CURRENT_SOURCE_DIR}/../rtexif ${EXTRA_INCDIR} ${GLIB2_INCLUDE_DIRS} ${GLIBMM_INCLUDE_DIRS}
    ${GTK_INCLUDE_DIRS} ${GTKMM_INCLUDE_DIRS} ${GIO_INCLUDE_DIRS} ${GIOMM_INCLUDE_DIRS} ${IPTCDATA_INCLUDE_DIRS} ${LCMS_INCLUDE_DIRS} ${GTHREAD_INCLUDE_DIRS} ${GOBJECT_INCLUDE_DIRS} )
link_directories (. ${CMAKE_CURRENT_SOURCE_DIR}/../rtexif ${EXTRA_LIBDIR} ${GLIB2_LIBRARY_DIRS} ${GLIBMM_LIBRARY_DIRS}
    ${GTK_LIBRARY_DIRS} ${GTKMM_LIBRARY_DIRS} ${GIO_LIBRARY_DIRS} ${GIOMM_LIBRARY_DIRS} ${IPTCDATA_LIBRARY_DIRS} ${LCMS_LIBRARY_DIRS} ${GTHREAD_LIBRARY_DIRS} ${GOBJECT_LIBRARY_DIRS})

    add_executable (rth ${BASESOURCEFILES})
ENDIF (WIN32)


set_target_properties (rth PROPERTIES COMPILE_FLAGS "-O3" OUTPUT_NAME rt)
target_link_libraries (rth rtengine ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${ZLIB_LIBRARIES} ${TIFF_LIBRARIES} ${EXTRA_LIB} ${GOBJECT_LIBRARIES} ${GTHREAD_LIBRARIES}
   ${GLIB2_LIBRARIES} ${GLIBMM_LIBRARIES} ${GTK_LIBRARIES} ${GTKMM_LIBRARIES} ${GIO_LIBRARIES} ${GIOMM_LIBRARIES} ${LCMS_LIBRARIES} ${IPTCDATA_LIBRARIES})

SET (INSTALL_SOURCE_PATH "${CMAKE_BINARY_DIR}/rtgui")

IF (WIN32)
    install (FILES ${INSTALL_SOURCE_PATH}/rth.exe DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../release
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
    install (FILES ${INSTALL_SOURCE_PATH}/../rawzor_win/rwz_sdk_s.dll DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../release
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
    install (FILES ${INSTALL_SOURCE_PATH}/../options.win DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../release
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ
	RENAME options)
ELSE (WIN32)
    install (FILES ${INSTALL_SOURCE_PATH}/rt DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../release
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
    install (FILES ../options.lin DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../release
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ
	RENAME options)
    install (FILES ../rtstart DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../release
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
    IF (CMAKE_SIZEOF_VOID_P EQUAL 4)
        install (FILES ../rawzor_lin32/librwz_sdk.so DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../release
            PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ
	        RENAME rwz_sdk.so)
    ELSEIF (CMAKE_SIZEOF_VOID_P EQUAL 8)
        install (FILES ../rawzor_lin64/librwz_sdk.so DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../release
            PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ
	        RENAME rwz_sdk.so)
    ENDIF (CMAKE_SIZEOF_VOID_P EQUAL 4)
ENDIF (WIN32)

