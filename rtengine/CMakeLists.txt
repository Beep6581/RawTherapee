
include_directories (${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../rtexif
    ${EXTRA_INCDIR} ${GTHREAD_INCLUDE_DIRS} ${GOBJECT_INCLUDE_DIRS} ${GLIB2_INCLUDE_DIRS}
    ${GLIBMM_INCLUDE_DIRS} ${IPTCDATA_INCLUDE_DIRS} ${LCMS_INCLUDE_DIRS})
link_directories (${CMAKE_CURRENT_SOURCE_DIR}/../rtexif ${EXTRA_LIBDIR} ${GTHREAD_LIBRARY_DIRS}
    ${GOBJECT_LIBRARY_DIRS} ${GLIB2_LIBRARY_DIRS} ${GLIBMM_LIBRARY_DIRS}
    ${IPTCDATA_LIBRARY_DIRS} ${LCMS_LIBRARY_DIRS})

set (RTENGINESOURCEFILES colortemp.cc curves.cc dcraw.cc gauss.cc iccstore.cc
    image8.cc image16.cc imagedata.cc imageio.cc improcfun.cc init.cc dcrop.cc
    loadinitial.cc procparams.cc rawimagesource.cc shmap.cc simpleprocess.cc refreshmap.cc
    stdimagesource.cc myfile.cc iccjpeg.c hlmultipliers.cc improccoordinator.cc
    helpers.cc processingjob.cc rtthumbnail.cc utils.cc bilateral2.cc)

if (BUILD_SHARED)
    add_library (rtengine SHARED ${RTENGINESOURCEFILES})
    # it may be nice to store library version too
    install (TARGETS rtengine DESTINATION ${LIBDIR})
else (BUILD_SHARED)
    add_library (rtengine STATIC ${RTENGINESOURCEFILES})
endif (BUILD_SHARED)

set_target_properties (rtengine PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -funroll-loops")

target_link_libraries (rtengine rtexif ${EXTRA_LIB} ${GOBJECT_LIBRARIES} ${GTHREAD_LIBRARIES}
    ${GLIB2_LIBRARIES} ${GLIBMM_LIBRARIES} ${LCMS_LIBRARIES} ${IPTCDATA_LIBRARIES}
    ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${TIFF_LIBRARIES} ${ZLIB_LIBRARIES})
