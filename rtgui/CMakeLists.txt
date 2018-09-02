# Common source files for both CLI and non-CLI execautables
set(CLISOURCEFILES
    alignedmalloc.cc
    custom-widgets/on-preview-editing/edit.cc
    main-cli.cc
    multilangmgr.cc
    options.cc
    paramsedited.cc
    pathutils.cc
    threadutils.cc
    )

set(NONCLISOURCEFILES
    alignedmalloc.cc
    batchqueue.cc
    batchqueuebuttonset.cc
    batchqueueentry.cc
    batchqueuepanel.cc
    batchtoolpanelcoord.cc
    bqentryupdater.cc
    browserfilter.cc
    cacheimagedata.cc
    cachemanager.cc
    cacorrection.cc
    clipboard.cc
    crophandler.cc
    cursormanager.cc
    custom-widgets/adjuster.cc
    custom-widgets/checkbox.cc
    custom-widgets/coloredbar.cc
    custom-widgets/coordinateadjuster.cc
    custom-widgets/curves/curveeditor.cc
    custom-widgets/curves/curveeditorgroup.cc
    custom-widgets/curves/diagonalcurveeditorsubgroup.cc
    custom-widgets/curves/flatcurveeditorsubgroup.cc
    custom-widgets/curves/mycurve.cc
    custom-widgets/curves/mydiagonalcurve.cc
    custom-widgets/curves/myflatcurve.cc
    custom-widgets/inspector.cc
    custom-widgets/lockablecolorpicker.cc
    custom-widgets/lwbutton.cc
    custom-widgets/lwbuttonset.cc
    custom-widgets/on-preview-editing/edit.cc
    custom-widgets/popupbutton.cc
    custom-widgets/popupcommon.cc
    custom-widgets/popuptogglebutton.cc
    custom-widgets/previewwindow.cc
    custom-widgets/profilestorecombobox.cc
    custom-widgets/rtimage.cc
    custom-widgets/shcselector.cc
    custom-widgets/splash.cc
    custom-widgets/thresholdadjuster.cc
    custom-widgets/thresholdselector.cc
    dirbrowser.cc
    editorpanel.cc
    eventmapper.cc
    exiffiltersettings.cc
    exportpanel.cc
    extprog.cc
    filebrowser.cc
    filebrowserentry.cc
    filecatalog.cc
    filepanel.cc
    filethumbnailbuttonset.cc
    filterpanel.cc
    guiutils.cc
    histogrampanel.cc
    history.cc
    iccprofilecreator.cc
    ilabel.cc
    imagearea.cc
    imageareapanel.cc
    indclippedpanel.cc
    labgrid.cc
    localcontrast.cc
    main.cc
    metadatapanel.cc
    multilangmgr.cc
    navigator.cc
    options.cc
    paramsedited.cc
    pathutils.cc
    placesbrowser.cc
    previewhandler.cc
    previewloader.cc
    previewmodepanel.cc
    profilepanel.cc
    recentbrowser.cc
    saveformatpanel.cc
    softlight.cc
    soundman.cc
    threadutils.cc
    thumbbrowserbase.cc
    thumbbrowserentrybase.cc
    thumbimageupdater.cc
    thumbnail.cc
    toolbar.cc
    toolpanel.cc
    toolpanelcoord.cc
    tools/bayerpreprocess.cc
    tools/bayerprocess.cc
    tools/bayerrawexposure.cc
    tools/blackwhite.cc
    tools/chmixer.cc
    tools/coarsepanel.cc
    tools/colorappearance.cc
    tools/colortoning.cc
    tools/crop.cc
    tools/cropwindow.cc
    tools/darkframe.cc
    tools/defringe.cc
    tools/dirpyrdenoise.cc
    tools/dirpyrequalizer.cc
    tools/distortion.cc
    tools/epd.cc
    tools/exifpanel.cc
    tools/fattaltonemap.cc
    tools/filmsimulation.cc
    tools/flatfield.cc
    tools/gradient.cc
    tools/hsvequalizer.cc
    tools/icmpanel.cc
    tools/impulsedenoise.cc
    tools/iptcpanel.cc
    tools/labcurve.cc
    tools/lensgeom.cc
    tools/lensprofile.cc
    tools/pcvignette.cc
    tools/perspective.cc
    tools/preprocess.cc
    tools/prsharpening.cc
    tools/rawcacorrection.cc
    tools/rawexposure.cc
    tools/resize.cc
    tools/retinex.cc
    tools/rgbcurves.cc
    tools/rotate.cc
    tools/sensorbayer.cc
    tools/sensorxtrans.cc
    tools/shadowshighlights.cc
    tools/sharpenedge.cc
    tools/sharpening.cc
    tools/sharpenmicro.cc
    tools/tonecurve.cc
    tools/vibrance.cc
    tools/vignetting.cc
    tools/wavelet.cc
    tools/whitebalance.cc
    tools/xtransprocess.cc
    tools/xtransrawexposure.cc
    window/dynamicprofilepanel.cc
    window/editwindow.cc
    window/partialpastedlg.cc
    window/preferences.cc
    window/renamedlg.cc
    window/rtwindow.cc
    window/saveasdlg.cc
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

if(OPENMP_FOUND AND NOT APPLE)
    set(EXTRA_LIB_RTGUI ${EXTRA_LIB_RTGUI} "atomic")
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
