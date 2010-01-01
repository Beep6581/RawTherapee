cmake_minimum_required(VERSION 2.6)
find_package(PkgConfig)

pkg_check_modules (GLIB2   glib-2.0>=2.16) 
pkg_check_modules (GLIBMM  glibmm-2.4>=2.16)
pkg_check_modules (SIGC    sigc++-2.0)

IF (WIN32)
    SET (EXTRA_LIBDIR "../lib; ../rawzor_win")
    SET (EXTRA_INCDIR "../winclude; ../rawzor_win")
    SET (EXTRA_LIB "ws2_32 ../lib/libiptcdata.a ../lib/libjpeg.a ../lib/libpng.a ../lib/libtiff.a ../lib/libz.a ../rawzor_win/rwz_sdk_s.a")
    ADD_DEFINITIONS (-DRAWZOR_SUPPORT)
ELSE (WIN32)
    IF (CMAKE_SIZEOF_VOID_P EQUAL 4)
        SET (EXTRA_INCDIR "../rawzor_lin32")
	SET (EXTRA_LIBDIR "../rawzor_lin32")
    ELSEIF (CMAKE_SIZEOF_VOID_P EQUAL 8)
        SET (EXTRA_INCDIR "../rawzor_lin64")
	SET (EXTRA_LIBDIR "../rawzor_lin64")
    ENDIF (CMAKE_SIZEOF_VOID_P EQUAL 4)
    SET (EXTRA_LIB "iptcdata -lrwz_sdk")
    ADD_DEFINITIONS (-DRAWZOR_SUPPORT)
ENDIF (WIN32)

include_directories (. ../rtexif ${EXTRA_INCDIR} ${GLIB2_INCLUDE_DIRS} ${GLIBMM_INCLUDE_DIRS})
link_directories (../rtexif ${EXTRA_LIBDIR} ${GLIB2_LIBRARY_DIRS} ${GLIBMM_LIBRARY_DIRS})

add_library (rtengine SHARED colortemp.cc curves.cc dcraw.cc gauss.cc iccstore.cc
    image8.cc image16.cc imagedata.cc imageio.cc improcfun.cc init.cc dcrop.cc
    loadinitial.cc procparams.cc rawimagesource.cc shmap.cc simpleprocess.cc refreshmap.cc
    stdimagesource.cc myfile.cc iccjpeg.c hlmultipliers.cc improccoordinator.cc
    helpers.cc processingjob.cc rtthumbnail.cc utils.cc hlmultipliers.cc bilateral2.cc)

set_target_properties (rtengine PROPERTIES COMPILE_FLAGS "-O3 -ffast-math -fexpensive-optimizations -funroll-loops -msse")

#add_executable (rtcmd rtetest.cc)
#add_executable (ppmap pparamsmap.cc)

target_link_libraries (rtengine rtexif lcms ${EXTRA_LIB} gthread-2.0 gobject-2.0 ${GLIB2_LIBRARIES} ${GLIBMM_LIBRARIES})
#target_link_libraries (rtcmd  rtengine)
#target_link_libraries (ppmap  rtengine)

IF (WIN32)
install (FILES librtengine.dll DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../release
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ         )
ELSE (WIN32)
install (FILES librtengine.so DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../release
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ         )
ENDIF (WIN32)
