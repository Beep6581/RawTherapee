cmake_minimum_required(VERSION 2.6)

# the default target is 'Debug'
set (CMAKE_BUILD_TYPE Debug CACHE STRING "One of: None Debug Release RelWithDebInfo MinSizeRel.")

# set the default target processor to 'native', i.e. the one of your own machine. The complete list of targets is in ProcessorTargets.cmake
# you can override the cpu list by setting PROC_TARGET_NUMBER to 0, but you have to provide at least
# the PROC_LABEL value (set to "undefined" by default) in order to provide downloadable builds,
# and evetually provide PROC_FLAGS if you want to manually set the processor specific flags
set (PROC_TARGET_NUMBER 2 CACHE STRING "Selected target processor from the list above (taken from ProcessorTargets.cmake)")

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

# adding the proc flags to the build flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PROC_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PROC_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${PROC_FLAGS}")

if (WIN32)
    set(CMAKE_PREFIX_PATH $ENV{MINGW_BASEPATH} $ENV{GTKMM_BASEPATH} CACHE STRING "Additional search paths")
endif (WIN32)

if (APPLE)
#    SET (CMAKE_OSX_ARCHITECTURES "i386;x86_64;" )
#    SET (CMAKE_TRY_COMPILE_OSX_ARCHITECTURES "i386;x86_64;" )
    SET (CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.5.sdk")
    SET (CMAKE_OSX_DEPLOYMENT_TARGET "10.5")
endif (APPLE)

option (BUILD_SHARED "Build rawtherapee with shared libraries" OFF)
option (WITH_RAWZOR "Build with Rawzor support" OFF)
option (WITH_MYFILE_MMAP "Build using memory mapped file" ON)
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
# NOTE: The new mechanism has been tested with BUILD_SHARED = OFF
if (WIN32)
    add_definitions (-DWIN32)
    add_definitions (-D_WIN32)
    if (MINGW)
        add_definitions (-D__MINGW32__)
    endif (MINGW)

    set (EXTRA_LIB "-lws2_32")
else (WIN32)
endif (WIN32)
pkg_check_modules (IPTCDATA REQUIRED libiptcdata)
pkg_check_modules (LCMS REQUIRED lcms<=1.99)
find_package (JPEG REQUIRED)
find_package (PNG REQUIRED)
find_package (TIFF REQUIRED)
find_package (ZLIB REQUIRED)

#generating AboutThisBuild.txt
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
    set (EXTRA_INCDIR ${EXTRA_INCDIR} "${CMAKE_CURRENT_SOURCE_DIR}/rawzor")
    add_definitions (-DRAWZOR_SUPPORT)
    if (WIN32)
        if (CMAKE_SIZEOF_VOID_P EQUAL 4)
            set (EXTRA_LIBDIR ${EXTRA_LIBDIR} "${CMAKE_CURRENT_SOURCE_DIR}/rawzor/win32")
            set (EXTRA_LIB ${EXTRA_LIB} "rwz_sdk.lib")
            install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/rawzor/win32/rwz_sdk.dll DESTINATION ${BINDIR}
                PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
        elseif (CMAKE_SIZEOF_VOID_P EQUAL 8)
            set (EXTRA_LIBDIR ${EXTRA_LIBDIR} "${CMAKE_CURRENT_SOURCE_DIR}/rawzor/win64")
            set (EXTRA_LIB ${EXTRA_LIB} "rwz_sdk.lib")
            install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/rawzor/win64/rwz_sdk.dll DESTINATION ${BINDIR}
                PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
        endif (CMAKE_SIZEOF_VOID_P EQUAL 4)
    elseif (APPLE)
        if (CMAKE_SIZEOF_VOID_P EQUAL 4)
            set (EXTRA_LIBDIR ${EXTRA_LIBDIR} "${CMAKE_CURRENT_SOURCE_DIR}/rawzor/mac32")
            set (EXTRA_LIB "-lrwz_sdk")
            install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/rawzor/mac32/rwz_sdk.dylib DESTINATION ${LIBDIR}
                PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
        elseif (CMAKE_SIZEOF_VOID_P EQUAL 8)
            set (EXTRA_LIBDIR ${EXTRA_LIBDIR} "${CMAKE_CURRENT_SOURCE_DIR}/rawzor/mac64")
            set (EXTRA_LIB "-lrwz_sdk_64")
            install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/rawzor/mac64/rwz_sdk_64.dylib DESTINATION ${LIBDIR}
                PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
        endif (CMAKE_SIZEOF_VOID_P EQUAL 4)
    else (WIN32)
        if (CMAKE_SIZEOF_VOID_P EQUAL 4)
            set (EXTRA_LIBDIR ${EXTRA_LIBDIR} "${CMAKE_CURRENT_SOURCE_DIR}/rawzor/lin32")
            set (EXTRA_LIB ${EXTRA_LIB} "-lrwz_sdk")
            install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/rawzor/lin32/librwz_sdk.so DESTINATION ${LIBDIR}
                PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ RENAME rwz_sdk.so)
        elseif (CMAKE_SIZEOF_VOID_P EQUAL 8)
            set (EXTRA_LIBDIR ${EXTRA_LIBDIR} "${CMAKE_CURRENT_SOURCE_DIR}/rawzor/lin64")
            set (EXTRA_LIB ${EXTRA_LIB} "-lrwz_sdk")
            install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/rawzor/lin64/librwz_sdk.so DESTINATION ${LIBDIR}
                PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ RENAME rwz_sdk.so)
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
