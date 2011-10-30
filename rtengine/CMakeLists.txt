
include_directories (${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../rtexif
    ${EXTRA_INCDIR} ${GTHREAD_INCLUDE_DIRS} ${GOBJECT_INCLUDE_DIRS} ${GLIB2_INCLUDE_DIRS}
    ${GLIBMM_INCLUDE_DIRS} ${IPTCDATA_INCLUDE_DIRS} ${LCMS_INCLUDE_DIRS}
	${CMAKE_CURRENT_SOURCE_DIR}/../rtgui ${GTKMM_INCLUDE_DIRS} ${GTK_INCLUDE_DIRS})
link_directories (${CMAKE_CURRENT_SOURCE_DIR}/../rtexif ${EXTRA_LIBDIR} ${GTHREAD_LIBRARY_DIRS}
    ${GOBJECT_LIBRARY_DIRS} ${GLIB2_LIBRARY_DIRS} ${GLIBMM_LIBRARY_DIRS}
    ${IPTCDATA_LIBRARY_DIRS} ${LCMS_LIBRARY_DIRS})

set (RTENGINESOURCEFILES safegtk.cc colortemp.cc curves.cc flatcurves.cc diagonalcurves.cc dcraw.cc iccstore.cc
    dfmanager.cc ffmanager.cc rawimage.cc image8.cc image16.cc imagefloat.cc imagedata.cc imageio.cc improcfun.cc init.cc dcrop.cc
    loadinitial.cc procparams.cc rawimagesource.cc demosaic_algos.cc shmap.cc simpleprocess.cc refreshmap.cc
    stdimagesource.cc myfile.cc iccjpeg.c hlmultipliers.cc improccoordinator.cc
    processingjob.cc rtthumbnail.cc utils.cc labimage.cc slicer.cc
    iplab2rgb.cc ipsharpen.cc iptransform.cc ipresize.cc ipvibrance.cc
	jpeg_memsrc.c jdatasrc.c
	PF_correct_RT.cc
    dirpyrLab_denoise.cc dirpyrLab_equalizer.cc dirpyr_equalizer.cc
    calc_distort.c
    klt/convolve.c klt/error.c klt/klt.c klt/klt_util.c klt/pnmio.c klt/pyramid.c klt/selectGoodFeatures.c
    klt/storeFeatures.c klt/trackFeatures.c klt/writeFeatures.c
    )

add_library (rtengine ${RTENGINESOURCEFILES})
#It may be nice to store library version too
IF (BUILD_SHARED_LIBS)
    install (TARGETS rtengine DESTINATION ${LIBDIR})
ENDIF (BUILD_SHARED_LIBS)

set_target_properties (rtengine PROPERTIES COMPILE_FLAGS "${RTENGINE_CXX_FLAGS}")

target_link_libraries (rtengine rtexif ${EXTRA_LIB} ${GOBJECT_LIBRARIES} ${GTHREAD_LIBRARIES}
    ${GLIB2_LIBRARIES} ${GLIBMM_LIBRARIES} ${LCMS_LIBRARIES} ${IPTCDATA_LIBRARIES}
    ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${TIFF_LIBRARIES} ${ZLIB_LIBRARIES})
