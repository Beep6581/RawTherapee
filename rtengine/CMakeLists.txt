
include_directories (${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../rtexif
    ${EXTRA_INCDIR} ${GTHREAD_INCLUDE_DIRS} ${GOBJECT_INCLUDE_DIRS} ${GLIB2_INCLUDE_DIRS}
    ${GLIBMM_INCLUDE_DIRS} ${IPTCDATA_INCLUDE_DIRS} ${LCMS_INCLUDE_DIRS})
link_directories (${CMAKE_CURRENT_SOURCE_DIR}/../rtexif ${EXTRA_LIBDIR} ${GTHREAD_LIBRARY_DIRS}
    ${GOBJECT_LIBRARY_DIRS} ${GLIB2_LIBRARY_DIRS} ${GLIBMM_LIBRARY_DIRS}
    ${IPTCDATA_LIBRARY_DIRS} ${LCMS_LIBRARY_DIRS})

set (RTENGINESOURCEFILES colortemp.cc curves.cc dcraw.cc filter.cc filterchain.cc filterchaingroup.cc
		filterfactory.cc iccjpeg.c iccstore.cc image16.cc image8.cc imagedata.cc imageio.cc imagesource.cc
		imageview.cc improcfuns.cc interactiveimprocimpl.cc matrix33.cc multiimage.cc myfile.cc procparams.cc
		rawimage.cc rawimagesource.cc settings.cc stdimagesource.cc processingjob.cc singleimproc.cc
		loadinitial.cc jdatasrc.c)
		
set (RTENGINEFILTERFILES filtchmixer.cc filtcoarse.cc filtcolorcurve.cc filtcolordenoise.cc filtdemosaic.cc
		filtinitial.cc filtlumacurve.cc filtlumadenoise.cc filtsh.cc filtsharpener.cc filttonecurve.cc
		filtwb.cc filthlrec.cc filtcsconv.cc filtresize.cc filttransform.cc)		

add_library (rtengine ${RTENGINESOURCEFILES} ${RTENGINEFILTERFILES})

#It may be nice to store library version too
IF (BUILD_SHARED_LIBS)
    install (TARGETS rtengine DESTINATION ${LIBDIR})
ENDIF (BUILD_SHARED_LIBS)

set_target_properties (rtengine PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -funroll-loops")

target_link_libraries (rtengine rtexif ${EXTRA_LIB} ${GOBJECT_LIBRARIES} ${GTHREAD_LIBRARIES}
    ${GLIB2_LIBRARIES} ${GLIBMM_LIBRARIES} ${LCMS_LIBRARIES} ${IPTCDATA_LIBRARIES}
    ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${TIFF_LIBRARIES} ${ZLIB_LIBRARIES})

add_executable (rtcmd rtcmd.cc)
target_link_libraries (rtcmd rtengine)
    