find_package(PkgConfig)

pkg_check_modules (GLIB2   REQUIRED glib-2.0>=2.16)
pkg_check_modules (GTHREAD REQUIRED gthread-2.0>=2.16)
pkg_check_modules (GOBJECT REQUIRED gobject-2.0>=2.16)
pkg_check_modules (GLIBMM  REQUIRED glibmm-2.4>=2.16)
pkg_check_modules (SIGC    REQUIRED sigc++-2.0)

IF (WIN32)
    SET (EXTRA_LIBDIR "${CMAKE_CURRENT_SOURCE_DIR}/../lib; ${CMAKE_CURRENT_SOURCE_DIR}/../rawzor_win")
    SET (EXTRA_INCDIR "${CMAKE_CURRENT_SOURCE_DIR}/../winclude; ${CMAKE_CURRENT_SOURCE_DIR}/../rawzor_win")
    SET (EXTRA_LIB "ws2_32 ${CMAKE_CURRENT_SOURCE_DIR}/../lib/libiptcdata.a ${CMAKE_CURRENT_SOURCE_DIR}/../lib/libjpeg.a ${CMAKE_CURRENT_SOURCE_DIR}/../lib/libpng.a ${CMAKE_CURRENT_SOURCE_DIR}/../lib/libtiff.a ${CMAKE_CURRENT_SOURCE_DIR}/../lib/libz.a")
    if (WITH_RAWZOR)
    ADD_DEFINITIONS (-DRAWZOR_SUPPORT)
    SET (EXTRA_LIB "${EXTRA_LIB} ${CMAKE_CURRENT_SOURCE_DIR}/../rawzor_win/rwz_sdk_s.a")
    endif (WITH_RAWZOR)
ELSE (WIN32)
    if (WITH_RAWZOR)
    IF (CMAKE_SIZEOF_VOID_P EQUAL 4)
        SET (EXTRA_INCDIR "${CMAKE_CURRENT_SOURCE_DIR}/../rawzor_lin32")
	SET (EXTRA_LIBDIR "${CMAKE_CURRENT_SOURCE_DIR}/../rawzor_lin32")
    ELSEIF (CMAKE_SIZEOF_VOID_P EQUAL 8)
        SET (EXTRA_INCDIR "${CMAKE_CURRENT_SOURCE_DIR}/../rawzor_lin64")
	SET (EXTRA_LIBDIR "${CMAKE_CURRENT_SOURCE_DIR}/../rawzor_lin64")
    ENDIF (CMAKE_SIZEOF_VOID_P EQUAL 4)
    SET (EXTRA_LIB "-lrwz_sdk")
    ADD_DEFINITIONS (-DRAWZOR_SUPPORT)
    endif (WITH_RAWZOR)
    pkg_check_modules (IPTCDATA REQUIRED libiptcdata)
    pkg_check_modules (LCMS REQUIRED lcms)
ENDIF (WIN32)

include_directories (${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../rtexif ${EXTRA_INCDIR} ${GTHREAD_INCLUDE_DIRS} ${GOBJECT_INCLUDE_DIRS} ${GLIB2_INCLUDE_DIRS} ${GLIBMM_INCLUDE_DIRS} ${IPTCDATA_INCLUDE_DIRS} ${LCMS_INCLUDE_DIRS})
link_directories (${CMAKE_CURRENT_SOURCE_DIR}/../rtexif ${EXTRA_LIBDIR} ${GTHREAD_LIBRARY_DIRS} ${GOBJECT_LIBRARY_DIRS} ${GLIB2_LIBRARY_DIRS} ${GLIBMM_LIBRARY_DIRS} ${IPTCDATA_LIBRARY_DIRS} ${LCMS_LIBRARY_DIRS})

add_library (rtengine SHARED colortemp.cc curves.cc dcraw.cc gauss.cc iccstore.cc
    image8.cc image16.cc imagedata.cc imageio.cc improcfun.cc init.cc dcrop.cc
    loadinitial.cc procparams.cc rawimagesource.cc shmap.cc simpleprocess.cc refreshmap.cc
    stdimagesource.cc myfile.cc iccjpeg.c hlmultipliers.cc improccoordinator.cc
    helpers.cc processingjob.cc rtthumbnail.cc utils.cc bilateral2.cc)

set_target_properties (rtengine PROPERTIES COMPILE_FLAGS "-O3 -ffast-math -fexpensive-optimizations -funroll-loops -msse")

#add_executable (rtcmd rtetest.cc)
#add_executable (ppmap pparamsmap.cc)

target_link_libraries (rtengine rtexif ${EXTRA_LIB} ${GOBJECT_LIBRARIES} ${GTHREAD_LIBRARIES} ${GLIB2_LIBRARIES} ${GLIBMM_LIBRARIES} ${LCMS_LIBRARIES} ${IPTCDATA_LIBRARIES})
#target_link_libraries (rtcmd  rtengine)
#target_link_libraries (ppmap  rtengine)

SET (INSTALL_SOURCE_PATH "${CMAKE_BINARY_DIR}/rtengine")


IF (WIN32)
install (FILES ${INSTALL_SOURCE_PATH}/librtengine.dll DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../release
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ         )
ELSE (WIN32)
install (FILES ${INSTALL_SOURCE_PATH}/librtengine.so DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../release
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ         )
ENDIF (WIN32)
