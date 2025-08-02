# Common source files for both CLI and non-CLI execautables
set(CLISOURCEFILES
    alignedmalloc.cc
    editcallbacks.cc
    main-cli.cc
    multilangmgr.cc
    options.cc
    paramsedited.cc
    pathutils.cc
    threadutils.cc
)

set(RT_GUI_TOOL_PANEL_SOURCES
    tools/bayerpreprocess.cc
    tools/bayerprocess.cc
    tools/bayerrawexposure.cc
    tools/blackwhite.cc
    tools/cacorrection.cc
    tools/chmixer.cc
    tools/colorappearance.cc
    tools/colortoning.cc
    tools/compressgamut.cc
    tools/crop.cc
    tools/darkframe.cc
    tools/defringe.cc
    tools/dehaze.cc
    tools/dirpyrdenoise.cc
    tools/dirpyrequalizer.cc
    tools/distortion.cc
    tools/epd.cc
    tools/fattaltonemap.cc
    tools/filmnegative.cc
    tools/filmsimulation.cc
    tools/flatfield.cc
    tools/framing.cc
    tools/gradient.cc
    tools/hsvequalizer.cc
    tools/icmpanel.cc
    tools/impulsedenoise.cc
    tools/labcurve.cc
    tools/lensgeom.cc
    tools/lensprofile.cc
    tools/localcontrast.cc
    tools/locallab.cc
    tools/locallabtools.cc
    tools/locallabtools2.cc
    tools/pcvignette.cc
    tools/pdsharpening.cc
    tools/perspective.cc
    tools/preprocess.cc
    tools/preprocesswb.cc
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
    tools/softlight.cc
    tools/spot.cc
    tools/tonecurve.cc
    tools/toneequalizer.cc
    tools/vibrance.cc
    tools/vignetting.cc
    tools/wavelet.cc
    tools/whitebalance.cc
    tools/xtransprocess.cc
    tools/xtransrawexposure.cc
)

set(NONCLISOURCEFILES
    adjuster.cc
    alignedmalloc.cc
    aspectratios.cc
    batchqueue.cc
    batchqueuebuttonset.cc
    batchqueueentry.cc
    batchqueuepanel.cc
    batchtoolpanelcoord.cc
    bqentryupdater.cc
    browserfilter.cc
    cacheimagedata.cc
    cachemanager.cc
    checkbox.cc
    clipboard.cc
    coarsepanel.cc
    coloredbar.cc
    colorpreview.cc
    controllines.cc
    controlspotpanel.cc
    coordinateadjuster.cc
    crophandler.cc
    cropwindow.cc
    cursormanager.cc
    curveeditor.cc
    curveeditorgroup.cc
    diagonalcurveeditorsubgroup.cc
    dirbrowser.cc
    dynamicprofilepanel.cc
    editbuffer.cc
    editcallbacks.cc
    editorpanel.cc
    editwidgets.cc
    editwindow.cc
    eventmapper.cc
    exiffiltersettings.cc
    exifpanel.cc
    exportpanel.cc
    externaleditorpreferences.cc
    extprog.cc
    filebrowser.cc
    filebrowserentry.cc
    filecatalog.cc
    filepanel.cc
    filethumbnailbuttonset.cc
    filterpanel.cc
    flatcurveeditorsubgroup.cc
    guiutils.cc
    hidpi.cc
    histogrampanel.cc
    history.cc
    iccprofilecreator.cc
    imagearea.cc
    imageareapanel.cc
    indclippedpanel.cc
    inspector.cc
    iptcpanel.cc
    labgrid.cc
    lockablecolorpicker.cc
    lwbutton.cc
    lwbuttonset.cc
    main.cc
    metadatapanel.cc
    multilangmgr.cc
    mycurve.cc
    mydiagonalcurve.cc
    myflatcurve.cc
    navigator.cc
    options.cc
    paramsedited.cc
    partialpastedlg.cc
    pathutils.cc
    placesbrowser.cc
    popupbutton.cc
    popupcommon.cc
    popuptogglebutton.cc
    preferences.cc
    previewhandler.cc
    previewloader.cc
    previewmodepanel.cc
    previewwindow.cc
    profilepanel.cc
    profilestorecombobox.cc
    recentbrowser.cc
    renamedlg.cc
    rtappchooserdialog.cc
    rtimage.cc
    rtscalable.cc
    rtsurface.cc
    rtwindow.cc
    saveasdlg.cc
    saveformatpanel.cc
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
    toolbar.cc
    toollocationpref.cc
    toolpanel.cc
    toolpanelcoord.cc
    zoompanel.cc
    ${RT_GUI_TOOL_PANEL_SOURCES}
)

if(APPLE)
    find_package(MacIntegration REQUIRED)
    # At the time of writing CMake has no module finder for gtkmacintegration so here we have it hard-coded, if installed via macports it should be in /opt/local/...
    set(EXTRA_LIB_RTGUI "${MacIntegration_LIBRARIES}")
    set(EXTRA_INCDIR "${EXTRA_INCDIR}" "${MacIntegration_INCLUDE_DIRS}")
endif()

if(WIN32)
    set(EXTRA_SRC_CLI myicon.rc)
    set(EXTRA_SRC_NONCLI myicon.rc)
    set(EXTRA_LIB_RTGUI winmm)
endif()

# Create new executables targets
add_executable(rth "${EXTRA_SRC_NONCLI}" "${NONCLISOURCEFILES}")
add_executable(rth-cli "${EXTRA_SRC_CLI}" "${CLISOURCEFILES}")

# Subdirectories that are not libraries themselves should be included
# like #include "your_subdirectory/header.h"
target_include_directories(rth PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(rth-cli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Add dependencies to executables targets
add_dependencies(rth UpdateInfo)
add_dependencies(rth-cli UpdateInfo)

# Create config.h which defines where data are stored
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/config.h")
# UpdateInfo and above creates new headers in build directory
target_include_directories(rth BEFORE PRIVATE ${CMAKE_BINARY_DIR}/rtgui)
target_include_directories(rth-cli BEFORE PRIVATE ${CMAKE_BINARY_DIR}/rtgui)

rt_setup_target(rth)
rt_setup_target(rth-cli)

#Define a target specific definition to use in code
target_compile_definitions(rth PUBLIC GUIVERSION)
target_compile_definitions(rth-cli PUBLIC CLIVERSION)

# Set executables targets properties, i.e. output filename and compile flags
# for "Debug" builds, open a console in all cases for Windows version
if((WIN32) AND NOT (CMAKE_BUILD_TYPE STREQUAL "Debug"))
    set_target_properties(rth PROPERTIES LINK_FLAGS "-mwindows")
endif()
set_target_properties(rth PROPERTIES OUTPUT_NAME rawtherapee)
set_target_properties(rth-cli PROPERTIES OUTPUT_NAME rawtherapee-cli)

target_include_directories(rth PRIVATE ${EXTRA_INCDIR})
target_include_directories(rth-cli PRIVATE ${EXTRA_INCDIR})

target_link_directories(rth PRIVATE ${EXTRA_LIBDIR})
target_link_directories(rth-cli PRIVATE ${EXTRA_LIBDIR})

# Add linked libraries dependencies to executables targets
target_link_libraries(rth PRIVATE
    rtengine
    PkgConfig::GTKMM
    ${EXTRA_LIB_RTGUI}
)
target_link_libraries(rth-cli PRIVATE
    rtengine
    PkgConfig::GTKMM
    ${EXTRA_LIB_RTGUI}
)

if(TCMALLOC_FOUND)
    target_link_libraries(rth PRIVATE TCMALLOC::TCMALLOC)
    target_link_libraries(rth-cli PRIVATE TCMALLOC::TCMALLOC)
endif()
if(USE_LIBCANBERRA)
    target_link_libraries(rth PRIVATE PkgConfig::CANBERRA)
endif()

# Install executables
install(TARGETS rth DESTINATION "${BINDIR}")
install(TARGETS rth-cli DESTINATION "${BINDIR}")
