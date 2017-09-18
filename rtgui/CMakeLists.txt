# Common source files for both CLI and non-CLI execautables
set(CLISOURCEFILES
    edit.cc
    main-cli.cc
    multilangmgr.cc
    options.cc
    paramsedited.cc
    pathutils.cc
    threadutils.cc
    )

set(NONCLISOURCEFILES
    adjuster.cc
    batchqueue.cc
    batchqueuebuttonset.cc
    batchqueueentry.cc
    batchqueuepanel.cc
    batchtoolpanelcoord.cc
    bayerpreprocess.cc
    bayerprocess.cc
    bayerrawexposure.cc
    blackwhite.cc
    bqentryupdater.cc
    browserfilter.cc
    cacheimagedata.cc
    cachemanager.cc
    cacorrection.cc
    checkbox.cc
    chmixer.cc
    clipboard.cc
    coarsepanel.cc
    colorappearance.cc
    coloredbar.cc
    colortoning.cc
    coordinateadjuster.cc
    crop.cc
    crophandler.cc
    cropwindow.cc
    cursormanager.cc
    curveeditor.cc
    curveeditorgroup.cc
    darkframe.cc
    defringe.cc
    diagonalcurveeditorsubgroup.cc
    dirbrowser.cc
    dirpyrdenoise.cc
    dirpyrequalizer.cc
    distortion.cc
    dynamicprofilepanel.cc
    edit.cc
    editorpanel.cc
    editwindow.cc
    epd.cc
    exiffiltersettings.cc
    exifpanel.cc
    exportpanel.cc
    extprog.cc
    filebrowser.cc
    filebrowserentry.cc
    filecatalog.cc
    filepanel.cc
    filethumbnailbuttonset.cc
    filmsimulation.cc
    filterpanel.cc
    flatcurveeditorsubgroup.cc
    flatfield.cc
    gradient.cc
    guiutils.cc
    histogrampanel.cc
    history.cc
    hsvequalizer.cc
    icmpanel.cc
    ilabel.cc
    imagearea.cc
    imageareapanel.cc
    impulsedenoise.cc
    indclippedpanel.cc
    inspector.cc
    iptcpanel.cc
    labcurve.cc
    lensgeom.cc
    lensprofile.cc
	locallab.cc
    lockablecolorpicker.cc
    lwbutton.cc
    lwbuttonset.cc
    main.cc
    multilangmgr.cc
    mycurve.cc
    mydiagonalcurve.cc
    myflatcurve.cc
    navigator.cc
    options.cc
    paramsedited.cc
    partialpastedlg.cc
    pathutils.cc
    pcvignette.cc
    perspective.cc
    placesbrowser.cc
    popupbutton.cc
    popupcommon.cc
    popuptogglebutton.cc
    preferences.cc
    preprocess.cc
    previewhandler.cc
    previewloader.cc
    previewmodepanel.cc
    previewwindow.cc
    profilepanel.cc
    profilestorecombobox.cc
    prsharpening.cc
    rawcacorrection.cc
    rawexposure.cc
    recentbrowser.cc
    renamedlg.cc
    resize.cc
    retinex.cc
    rgbcurves.cc
    rotate.cc
    rtimage.cc
    rtwindow.cc
    saveasdlg.cc
    saveformatpanel.cc
    sensorbayer.cc
    sensorxtrans.cc
    shadowshighlights.cc
    sharpenedge.cc
    sharpening.cc
    sharpenmicro.cc
    shcselector.cc
    soundman.cc
    splash.cc
    threadutils.cc
    thresholdadjuster.cc
    thresholdselector.cc
    thumbbrowserbase.cc
    thumbbrowserentrybase.cc
    thumbimageupdater.cc
    thumbnail.cc
    tonecurve.cc
    toolbar.cc
    toolpanel.cc
    toolpanelcoord.cc
    vibrance.cc
    vignetting.cc
    wavelet.cc
    whitebalance.cc
    xtransprocess.cc
    xtransrawexposure.cc
    zoompanel.cc
    )

include_directories(BEFORE "${CMAKE_CURRENT_BINARY_DIR}")

if(APPLE)
    find_package(MacIntegration REQUIRED)
    # At the time of writing CMake has no module finder for gtkmacintegration so here we have it hard-coded, if installed via macports it should be in /opt/local/...
    set(EXTRA_LIB_RTGUI ${MacIntegration_LIBRARIES})
    set(EXTRA_INCDIR ${EXTRA_INCDIR} ${MacIntegration_INCLUDE_DIRS})
endif()

if(WIN32)
    set(EXTRA_SRC_CLI myicon.rc)
    set(EXTRA_SRC_NONCLI myicon.rc windirmonitor.cc)
    set(EXTRA_LIB_RTGUI winmm)
    include_directories(${EXTRA_INCDIR}
        ${GIOMM_INCLUDE_DIRS}
        ${GIO_INCLUDE_DIRS}
        ${GLIB2_INCLUDE_DIRS}
        ${GLIBMM_INCLUDE_DIRS}
        ${GTKMM_INCLUDE_DIRS}
        ${GTK_INCLUDE_DIRS}
        ${LENSFUN_INCLUDE_DIRS}
        )
    link_directories(. "${PROJECT_SOURCE_DIR}/rtexif"
        ${EXTRA_LIBDIR}
        ${GIOMM_LIBRARY_DIRS}
        ${GIO_LIBRARY_DIRS}
        ${GLIB2_LIBRARY_DIRS}
        ${GLIBMM_LIBRARY_DIRS}
        ${GTKMM_LIBRARY_DIRS}
        ${GTK_LIBRARY_DIRS}
        ${LENSFUN_LIBRARY_DIRS}
        )
else()
    include_directories(${EXTRA_INCDIR}
        ${CANBERRA-GTK_INCLUDE_DIRS}
        ${EXPAT_INCLUDE_DIRS}
        ${FFTW3F_LIBRARY_DIRS}
        ${GIOMM_INCLUDE_DIRS}
        ${GIO_INCLUDE_DIRS}
        ${GLIB2_INCLUDE_DIRS}
        ${GLIBMM_INCLUDE_DIRS}
        ${GOBJECT_INCLUDE_DIRS}
        ${GTHREAD_INCLUDE_DIRS}
        ${GTKMM_INCLUDE_DIRS}
        ${GTK_INCLUDE_DIRS}
        ${IPTCDATA_INCLUDE_DIRS}
        ${LCMS_INCLUDE_DIRS}
        ${LENSFUN_INCLUDE_DIRS}
        )
    link_directories(${EXTRA_LIBDIR}
        ${CANBERRA-GTK_LIBRARY_DIRS}
        ${EXPAT_LIBRARY_DIRS}
        ${FFTW3F_LIBRARY_DIRS}
        ${GIOMM_LIBRARY_DIRS}
        ${GIO_LIBRARY_DIRS}
        ${GLIB2_LIBRARY_DIRS}
        ${GLIBMM_LIBRARY_DIRS}
        ${GOBJECT_LIBRARY_DIRS}
        ${GTHREAD_LIBRARY_DIRS}
        ${GTKMM_LIBRARY_DIRS}
        ${GTK_LIBRARY_DIRS}
        ${IPTCDATA_LIBRARY_DIRS}
        ${LCMS_LIBRARY_DIRS}
        ${LENSFUN_LIBRARY_DIRS}
        )
endif()

# Create config.h which defines where data are stored
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/config.h")

# Create new executables targets
add_executable(rth ${EXTRA_SRC_NONCLI} ${NONCLISOURCEFILES})
add_executable(rth-cli ${EXTRA_SRC_CLI} ${CLISOURCEFILES})

# Add dependencies to executables targets
add_dependencies(rth UpdateInfo)
add_dependencies(rth-cli UpdateInfo)

# Set executables targets properties, i.e. output filename and compile flags
# for "Debug" builds, open a console in all cases for Windows version
if((WIN32) AND NOT(UPPER_CMAKE_BUILD_TYPE STREQUAL "DEBUG"))
    set_target_properties(rth PROPERTIES LINK_FLAGS "-mwindows")
endif()
set_target_properties(rth PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS}" OUTPUT_NAME rawtherapee)
set_target_properties(rth-cli PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS}" OUTPUT_NAME rawtherapee-cli)

# Add linked libraries dependencies to executables targets
target_link_libraries(rth rtengine
    ${CANBERRA-GTK_LIBRARIES}
    ${EXPAT_LIBRARIES}
    ${EXTRA_LIB_RTGUI}
    ${FFTW3F_LIBRARIES}
    ${GIOMM_LIBRARIES}
    ${GIO_LIBRARIES}
    ${GLIB2_LIBRARIES}
    ${GLIBMM_LIBRARIES}
    ${GOBJECT_LIBRARIES}
    ${GTHREAD_LIBRARIES}
    ${GTKMM_LIBRARIES}
    ${GTK_LIBRARIES}
    ${IPTCDATA_LIBRARIES}
    ${JPEG_LIBRARIES}
    ${LCMS_LIBRARIES}
    ${PNG_LIBRARIES}
    ${TIFF_LIBRARIES}
    ${ZLIB_LIBRARIES}
    ${LENSFUN_LIBRARIES}
    )

target_link_libraries(rth-cli rtengine
    ${CAIROMM_LIBRARIES}
    ${EXPAT_LIBRARIES}
    ${EXTRA_LIB_RTGUI}
    ${FFTW3F_LIBRARIES}
    ${GIOMM_LIBRARIES}
    ${GIO_LIBRARIES}
    ${GLIB2_LIBRARIES}
    ${GLIBMM_LIBRARIES}
    ${GOBJECT_LIBRARIES}
    ${GTHREAD_LIBRARIES}
    ${IPTCDATA_LIBRARIES}
    ${JPEG_LIBRARIES}
    ${LCMS_LIBRARIES}
    ${PNG_LIBRARIES}
    ${TIFF_LIBRARIES}
    ${ZLIB_LIBRARIES}
    ${LENSFUN_LIBRARIES}
    )

# Install executables
install(TARGETS rth DESTINATION ${BINDIR})
install(TARGETS rth-cli DESTINATION ${BINDIR})
