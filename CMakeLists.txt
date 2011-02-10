cmake_minimum_required(VERSION 2.6)

# the default target is 'Debug'
set (CMAKE_BUILD_TYPE Debug CACHE STRING "One of: None Debug Release RelWithDebInfo MinSizeRel.")

# By default, we don't use specific processor target, so  PROC_TARGET_NUMBER is set to 0. If can specify other values to select specific
# processor targets, which list can be found in ProcessorTargets.cmake.
set (PROC_TARGET_NUMBER 0 CACHE STRING "Selected target processor from the list above (taken from ProcessorTargets.cmake)")

# The following line set special compilation flags for RTEngine, and will be added to CMAKE_CXX_FLAGS
# It were moved away from rtengine/CMakefiles.txt, because some users may want to remove -ffast_math :
# this flag speeds up the floating-point operations, but with a little bite less precisions. This default value
# gives the same result/behaviour as before.
set (RTENGINE_CXX_FLAGS "-ffast-math -funroll-loops" CACHE STRING "Special compilation flags for RTEngine")

#loading the processor targets list
set (PROC_LABEL "undefined")
set (PROC_FLAGS "")
if (PROC_TARGET_NUMBER GREATER 0)
	include (ProcessorTargets.cmake)
	set (PROC_LABEL ${PROC_TARGET_${PROC_TARGET_NUMBER}_LABEL})
	set (PROC_FLAGS ${PROC_TARGET_${PROC_TARGET_NUMBER}_FLAGS})
endif (PROC_TARGET_NUMBER GREATER 0)

# if it exists, the PROC_FORCED_LABEL value is loaded in PROC_LABEL to override the one of ProcessorTargets
if (DEFINED PROC_FORCED_LABEL)
  set (PROC_LABEL ${PROC_FORCED_LABEL})
endif (DEFINED PROC_FORCED_LABEL)

# adding the proc flags to the build flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PROC_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PROC_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${PROC_FLAGS}")

if (APPLE)
#    SET (CMAKE_OSX_ARCHITECTURES "i386;x86_64;" )
#    SET (CMAKE_TRY_COMPILE_OSX_ARCHITECTURES "i386;x86_64;" )
    SET (CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.5.sdk")
    SET (CMAKE_OSX_DEPLOYMENT_TARGET "10.5")
endif (APPLE)

option (BUILD_SHARED "Build rawtherapee with shared libraries" OFF)
option (WITH_RAWZOR "Build with Rawzor support" OFF)
option (WITH_MYFILE_MMAP "Build using memory mapped file" OFF)
option (OPTION_OMP "Build with OpenMP support" ON)

# set install directories
if (NOT DEFINED DATADIR)
    if (WIN32 OR APPLE)
        set (DATADIR ${CMAKE_CURRENT_SOURCE_DIR}/release)
    else (WIN32 OR APPLE)
        set (DATADIR ${CMAKE_INSTALL_PREFIX}/share/rawtherapee)
    endif (WIN32 OR APPLE)
endif (NOT DEFINED DATADIR)

if (NOT DEFINED BINDIR)
    if (WIN32 OR APPLE)
        set (BINDIR  ${CMAKE_CURRENT_SOURCE_DIR}/release)
    else (WIN32 OR APPLE)
        set (BINDIR  ${CMAKE_INSTALL_PREFIX}/bin)
    endif (WIN32 OR APPLE)
endif (NOT DEFINED BINDIR)

if (NOT DEFINED LIBDIR)
    if (WIN32 OR APPLE)
        set (LIBDIR ${CMAKE_CURRENT_SOURCE_DIR}/release)
    else (WIN32 OR APPLE)
        # Respect CMAKE_INSTALL_LIBDIR if set
        if (DEFINED CMAKE_INSTALL_LIBDIR)
            set (LIBDIR ${CMAKE_INSTALL_LIBDIR})
        else (DEFINED CMAKE_INSTALL_LIBDIR)
            set (LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
        endif (DEFINED CMAKE_INSTALL_LIBDIR)
    endif (WIN32 OR APPLE)
endif (NOT DEFINED LIBDIR)

if (NOT DEFINED DOCDIR)
    if (WIN32 OR APPLE)
        set (DOCDIR ${CMAKE_CURRENT_SOURCE_DIR}/release/doc)
    else (WIN32 OR APPLE)
        set (DOCDIR ${CMAKE_INSTALL_PREFIX}/share/doc)
    endif (WIN32 OR APPLE)
endif (NOT DEFINED DOCDIR)

if (NOT DEFINED CREDITSDIR)
    if (WIN32 OR APPLE)
        set (CREDITSDIR ${CMAKE_CURRENT_SOURCE_DIR}/release)
    else (WIN32 OR APPLE)
        set (CREDITSDIR ${CMAKE_INSTALL_PREFIX}/share/doc)
    endif (WIN32 OR APPLE)
endif (NOT DEFINED CREDITSDIR)

if (NOT DEFINED LICENCEDIR)
    if (WIN32 OR APPLE)
        set (LICENCEDIR ${CMAKE_CURRENT_SOURCE_DIR}/release)
    else (WIN32 OR APPLE)
        set (LICENCEDIR ${CMAKE_INSTALL_PREFIX}/share/doc)
    endif (WIN32 OR APPLE)
endif (NOT DEFINED LICENCEDIR)

# check for libraries
find_package(PkgConfig)
pkg_check_modules (GTK     REQUIRED gtk+-2.0>=2.12)
pkg_check_modules (GLIB2   REQUIRED glib-2.0>=2.16)
pkg_check_modules (GLIBMM  REQUIRED glibmm-2.4>=2.16)
pkg_check_modules (GTKMM   REQUIRED gtkmm-2.4>=2.12)
pkg_check_modules (GIO     REQUIRED gio-2.0>=2.16)
pkg_check_modules (GIOMM   REQUIRED giomm-2.4>=2.12)
pkg_check_modules (GTHREAD REQUIRED gthread-2.0>=2.16)
pkg_check_modules (GOBJECT REQUIRED gobject-2.0>=2.16)
pkg_check_modules (SIGC    REQUIRED sigc++-2.0)
# NOTE: dependencies should be handled by pkg_check_modules and FIND_PACKAGE
# on windows too but I don't want to break current build chain
if (WIN32)
    set (EXTRA_LIBDIR "${CMAKE_CURRENT_SOURCE_DIR}/lib")
    set (EXTRA_INCDIR "${CMAKE_CURRENT_SOURCE_DIR}/winclude")
    set (EXTRA_LIB "-lws2_32")
    set (IPTCDATA_LIBRARIES iptcdata)
    set (LCMS_LIBRARIES liblcms.a)
    set (JPEG_LIBRARIES libjpeg.a)
    set (PNG_LIBRARIES libpng.a)
    set (TIFF_LIBRARIES libtiff.a)
    set (ZLIB_LIBRARIES libz.a)

    add_definitions (-DWIN32)
    add_definitions (-D_WIN32)
    if (MINGW)
        add_definitions (-D__MINGW32__)
    endif (MINGW)
else (WIN32)
    pkg_check_modules (IPTCDATA REQUIRED libiptcdata)
    pkg_check_modules (LCMS    REQUIRED lcms)
    find_package (JPEG REQUIRED)
    find_package (PNG REQUIRED)
    find_package (TIFF REQUIRED)
    find_package (ZLIB REQUIRED)
endif (WIN32)

# set the bit number information of the platform
if (CMAKE_SIZEOF_VOID_P EQUAL 4)
    set(PROC_BIT_DEPTH 32 bits)
elseif (CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(PROC_BIT_DEPTH 64 bits)
endif (CMAKE_SIZEOF_VOID_P EQUAL 4)

#generating AboutThisBuild.txt
set (BUILD_TYPE CMAKE_BUILD_TYPE)
if (WIN32)
    include (About-Windows.cmake)
elseif (APPLE)
    include (About-Apple.cmake)
else (WIN32)
    include (About-Linux.cmake)
endif (WIN32)
add_dependencies(AboutFile Debug Release MinSizeRel RelWithDebInfo)

# link rawzor
if (WITH_RAWZOR)
    if (WIN32)
        set (EXTRA_INCDIR ${EXTRA_INCDIR} "${CMAKE_CURRENT_SOURCE_DIR}/rawzor_win")
        set (EXTRA_LIBDIR ${EXTRA_LIBDIR} "${CMAKE_CURRENT_SOURCE_DIR}/rawzor_win")
        set (EXTRA_LIB ${EXTRA_LIB} "${CMAKE_CURRENT_SOURCE_DIR}/rawzor_win/rwz_sdk_s.a")
        add_definitions (-DRAWZOR_SUPPORT)
        install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/rawzor_win/rwz_sdk_s.dll DESTINATION ${BINDIR}
            PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
    elseif (APPLE)
        set (EXTRA_INCDIR "${CMAKE_CURRENT_SOURCE_DIR}/rawzor_mac")
        set (EXTRA_LIBDIR "${CMAKE_CURRENT_SOURCE_DIR}/rawzor_mac")
        set (EXTRA_LIB "-lrwz_sdk_64")
        add_definitions (-DRAWZOR_SUPPORT)
        install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/rawzor_mac/librwz_sdk_64.dylib DESTINATION ${BINDIR}
            PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ RENAME rwz_sdk_64.dylib)
        if (BUILD_SHARED)
            install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/rawzor_mac/librwz_sdk_64.dylib DESTINATION ${LIBDIR}
                PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ RENAME rwz_sdk_64.dylib)
        endif (BUILD_SHARED)
    else (WIN32)
        if (CMAKE_SIZEOF_VOID_P EQUAL 4)
            set (EXTRA_INCDIR ${EXTRA_INCDIR} "${CMAKE_CURRENT_SOURCE_DIR}/rawzor_lin32")
            set (EXTRA_LIBDIR ${EXTRA_LIBDIR} "${CMAKE_CURRENT_SOURCE_DIR}/rawzor_lin32")
            set (EXTRA_LIB ${EXTRA_LIB} "-lrwz_sdk")
            add_definitions (-DRAWZOR_SUPPORT)
            install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/rawzor_lin32/librwz_sdk.so DESTINATION ${BINDIR}
                PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ RENAME rwz_sdk.so)
            if (BUILD_SHARED)
                install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/rawzor_lin32/librwz_sdk.so DESTINATION ${LIBDIR}
                    PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ RENAME rwz_sdk.so)
            endif (BUILD_SHARED)
        elseif (CMAKE_SIZEOF_VOID_P EQUAL 8)
            set (EXTRA_INCDIR ${EXTRA_INCDIR} "${CMAKE_CURRENT_SOURCE_DIR}/rawzor_lin64")
            set (EXTRA_LIBDIR ${EXTRA_LIBDIR} "${CMAKE_CURRENT_SOURCE_DIR}/rawzor_lin64")
            set (EXTRA_LIB ${EXTRA_LIB} "-lrwz_sdk")
            add_definitions (-DRAWZOR_SUPPORT)
            install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/rawzor_lin64/librwz_sdk.so DESTINATION ${BINDIR}
                PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ RENAME rwz_sdk.so)
            if (BUILD_SHARED)
                install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/rawzor_lin64/librwz_sdk.so DESTINATION ${LIBDIR}
                    PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ RENAME rwz_sdk.so)
            endif (BUILD_SHARED)
        endif (CMAKE_SIZEOF_VOID_P EQUAL 4)
    endif (WIN32)
endif (WITH_RAWZOR)

if (WITH_MYFILE_MMAP)
	add_definitions (-DMYFILE_MMAP)
endif (WITH_MYFILE_MMAP)

if (OPTION_OMP)
	find_package(OpenMP)
	if (OPENMP_FOUND)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif (OPENMP_FOUND)
endif (OPTION_OMP)

if (UNIX)
    install (PROGRAMS rtstart DESTINATION ${BINDIR})
endif (UNIX)

install (FILES AUTHORS.txt DESTINATION ${CREDITSDIR})
install (FILES LICENSE.txt DESTINATION ${LICENCEDIR})
install (FILES AboutThisBuild.txt DESTINATION ${BINDIR})

add_subdirectory (rtexif)
add_subdirectory (rtengine)
add_subdirectory (rtgui)
add_subdirectory (rtdata)
