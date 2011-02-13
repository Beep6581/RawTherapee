
include_directories (${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../rtexif
    ${EXTRA_INCDIR} ${LCMS_INCLUDE_DIRS} ${EXIV2_INCLUDE_DIRS}
    ${LIBRAW_INCLUDE_DIRS} ${FREEIMAGE_INCLUDE_DIRS} ${TOOLKIT_INCLUDE_DIRS})
link_directories (${CMAKE_CURRENT_SOURCE_DIR}/../rtexif ${EXTRA_LIBDIR} 
	${LCMS_LIBRARY_DIRS} ${EXIV2_LIBRARY_DIRS} ${LIBRAW_LIBRARY_DIRS} ${TOOLKIT_LINK_DIRS})

set (RTENGINESOURCEFILES  image.cc multiimage.cc procparams.cc colortemp.cc curves.cc filter.cc filterchain.cc filterchaingroup.cc
		filterfactory.cc iccstore.cc imagedata.cc imagesource.cc
		imageview.cc improcfuns.cc interactiveimprocimpl.cc matrix33.cc  
		rawimage.cc rawimagesource.cc settings.cc stdimagesource.cc processingjob.cc singleimproc.cc
		loadinitial.cc bilateral.cc)
		
set (RTENGINEFILTERFILES filtchmixer.cc filtcoarse.cc filtcolorcurve.cc filtcolordenoise.cc filtdemosaic.cc
		filtinitial.cc filtlumacurve.cc filtlumadenoise.cc filtsh.cc filtsharpener.cc filttonecurve.cc
		filtwb.cc filthlrec.cc filtcsconv.cc filtresize.cc filttransform.cc)		

add_library (rtengine ${RTENGINESOURCEFILES} ${RTENGINEFILTERFILES})

#It may be nice to store library version too
IF (BUILD_SHARED_LIBS)
    install (TARGETS rtengine DESTINATION ${LIBDIR})
ENDIF (BUILD_SHARED_LIBS)

target_link_libraries (rtengine ${EXTRA_LIB} ${LCMS_LIBRARIES}  ${EXIV2_LIBRARIES}
    ${LIBRAW_LIBRARIES} ${FREEIMAGE_LIBRARIES} ${TOOLKIT_LINK_LIBS})

add_executable (rtcmd rtcmd.cc)
target_link_libraries (rtcmd rtengine)

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)
    
