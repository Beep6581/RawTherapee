
set (BASESOURCEFILES
    editwindow.cc batchtoolpanelcoord.cc paramsedited.cc cropwindow.cc previewhandler.cc previewwindow.cc navigator.cc indclippedpanel.cc filterpanel.cc
    cursormanager.cc rtwindow.cc renamedlg.cc recentbrowser.cc placesbrowser.cc filepanel.cc editorpanel.cc batchqueuepanel.cc
    ilabel.cc thumbbrowserbase.cc adjuster.cc filebrowserentry.cc filebrowser.cc filethumbnailbuttonset.cc
    cachemanager.cc cacheimagedata.cc shcselector.cc perspective.cc
    clipboard.cc thumbimageupdater.cc bqentryupdater.cc lensgeom.cc
    coarsepanel.cc cacorrection.cc colorshift.cc hlrec.cc chmixer.cc
    colorboost.cc resize.cc icmpanel.cc crop.cc shadowshighlights.cc
    colordenoise.cc
    impulsedenoise.cc dirpyrdenoise.cc
    exifpanel.cc toolpanel.cc
    sharpening.cc
    whitebalance.cc   vignetting.cc rotate.cc distortion.cc
    crophandler.cc dirbrowser.cc
    curveeditor.cc curveeditorgroup.cc diagonalcurveeditorsubgroup.cc flatcurveeditorsubgroup.cc
    filecatalog.cc
    previewloader.cc
    histogrampanel.cc history.cc  imagearea.cc
    imageareapanel.cc iptcpanel.cc labcurve.cc lumadenoise.cc main.cc
    multilangmgr.cc mycurve.cc myflatcurve.cc mydiagonalcurve.cc options.cc
    preferences.cc profilepanel.cc safegtk.cc saveasdlg.cc
    saveformatpanel.cc splash.cc
    thumbnail.cc tonecurve.cc toolbar.cc
    guiutils.cc zoompanel.cc toolpanelcoord.cc
    thumbbrowserentrybase.cc batchqueueentry.cc
    batchqueue.cc lwbutton.cc lwbuttonset.cc
    batchqueuebuttonset.cc browserfilter.cc exiffiltersettings.cc
    profilestore.cc partialpastedlg.cc rawprocess.cc preprocess.cc
    darkframe.cc rawcacorrection.cc rawexposure.cc
    equalizer.cc dirpyrequalizer.cc hsvequalizer.cc defringe.cc
    popupcommon.cc popupbutton.cc popuptogglebutton.cc)

if (WIN32)
    set (EXTRA_SRC windirmonitor.cc myicon.rc)
    include_directories ( ../rtengine ${CMAKE_CURRENT_BINARY_DIR} . ../rtexif ${EXTRA_INCDIR} ${GLIB2_INCLUDE_DIRS} ${GLIBMM_INCLUDE_DIRS}
        ${GTK_INCLUDE_DIRS} ${GTKMM_INCLUDE_DIRS} ${GIO_INCLUDE_DIRS} ${GIOMM_INCLUDE_DIRS})
    link_directories (. ../rtexif ${EXTRA_LIBDIR} ${GLIB2_LIBRARY_DIRS} ${GLIBMM_LIBRARY_DIRS}
        ${GTK_LIBRARY_DIRS} ${GTKMM_LIBRARY_DIRS} ${GIO_LIBRARY_DIRS} ${GIOMM_LIBRARY_DIRS})
    #set_target_properties (rth PROPERTIES LINK_FLAGS "-mwindows")
else (WIN32)
    include_directories (${CMAKE_CURRENT_SOURCE_DIR}/../rtengine ${CMAKE_CURRENT_BINARY_DIR} . ${CMAKE_CURRENT_SOURCE_DIR}/../rtexif ${EXTRA_INCDIR} ${GLIB2_INCLUDE_DIRS} ${GLIBMM_INCLUDE_DIRS}
        ${GTK_INCLUDE_DIRS} ${GTKMM_INCLUDE_DIRS} ${GIO_INCLUDE_DIRS} ${GIOMM_INCLUDE_DIRS} ${IPTCDATA_INCLUDE_DIRS} ${LCMS_INCLUDE_DIRS} ${GTHREAD_INCLUDE_DIRS} ${GOBJECT_INCLUDE_DIRS} )
    link_directories (${EXTRA_LIBDIR} ${GLIB2_LIBRARY_DIRS} ${GLIBMM_LIBRARY_DIRS}
        ${GTK_LIBRARY_DIRS} ${GTKMM_LIBRARY_DIRS} ${GIO_LIBRARY_DIRS} ${GIOMM_LIBRARY_DIRS} ${IPTCDATA_LIBRARY_DIRS} ${LCMS_LIBRARY_DIRS} ${GTHREAD_LIBRARY_DIRS} ${GOBJECT_LIBRARY_DIRS})
    # create config.h which defines where data are stored
    configure_file (${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)
endif (WIN32)

add_executable (rth ${EXTRA_SRC} ${BASESOURCEFILES})

set_target_properties (rth PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS}" OUTPUT_NAME rt)
#target_link_libraries (rth rtengine ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${ZLIB_LIBRARIES} ${TIFF_LIBRARIES} ${EXTRA_LIB} ${GOBJECT_LIBRARIES} ${GTHREAD_LIBRARIES}
#   ${GLIB2_LIBRARIES} ${GLIBMM_LIBRARIES} ${GTK_LIBRARIES} ${GTKMM_LIBRARIES} ${GIO_LIBRARIES} ${GIOMM_LIBRARIES} ${LCMS_LIBRARIES} ${IPTCDATA_LIBRARIES})
target_link_libraries (rth rtengine ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${ZLIB_LIBRARIES} ${TIFF_LIBRARIES} ${GOBJECT_LIBRARIES} ${GTHREAD_LIBRARIES}
   ${GLIB2_LIBRARIES} ${GLIBMM_LIBRARIES} ${GTK_LIBRARIES} ${GTKMM_LIBRARIES} ${GIO_LIBRARIES} ${GIOMM_LIBRARIES} ${LCMS_LIBRARIES} ${IPTCDATA_LIBRARIES})
install (TARGETS rth DESTINATION ${BINDIR})

