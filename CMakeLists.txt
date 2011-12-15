if (WIN32)
  cmake_minimum_required(VERSION 2.8.4)
  cmake_policy(SET CMP0015 OLD)
else (WIN32)
  cmake_minimum_required(VERSION 2.6)
endif (WIN32)

# the default target is 'Debug'
if (CMAKE_BUILD_TYPE STREQUAL "")
  set (CMAKE_BUILD_TYPE Debug CACHE STRING "One of: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif ()

string (TOUPPER ${CMAKE_BUILD_TYPE} UPPER_CMAKE_BUILD_TYPE)

if (UPPER_CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  add_definitions (-D_DEBUG)
endif ()

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

set (CACHE_NAME_SUFFIX "" CACHE STRING "RawTherapee's cache folder suffix (leave empty to use the default suffix, i.e. latesttag)")

# By default, we don't use specific processor target, so  PROC_TARGET_NUMBER is set to 0. If can specify other values to select specific
# processor targets, which list can be found in ProcessorTargets.cmake.
set (PROC_TARGET_NUMBER 0 CACHE STRING "Selected target processor from the list above (taken from ProcessorTargets.cmake)")

# The following line set special compilation flags for RTEngine, and will be added to CMAKE_CXX_FLAGS
# It were moved away from rtengine/CMakefiles.txt, because some users may want to remove -ffast_math :
# this flag speeds up the floating-point operations, but with a little bite less precisions. This default value
# gives the same result/behaviour as before.
set (RTENGINE_CXX_FLAGS "-ffast-math -funroll-loops" CACHE STRING "Special compilation flags for RTEngine")

#loading the processor targets list
include (ProcessorTargets.cmake)
set (PROC_LABEL "undefined" CACHE STRING "Target processor label, unused if PROC_TARGET_NUMBER = 0 or 2")
set (PROC_FLAGS "" CACHE STRING "Target processor related build/link flags")
if ((NOT(PROC_TARGET_NUMBER NOT EQUAL 0)) AND (NOT(PROC_TARGET_NUMBER EQUAL 2)))
	set (PROC_LABEL ${PROC_TARGET_${PROC_TARGET_NUMBER}_LABEL})
	set (PROC_FLAGS ${PROC_TARGET_${PROC_TARGET_NUMBER}_FLAGS})
endif ()

# if it exists, the PROC_FORCED_LABEL value is loaded in PROC_LABEL to override the one of ProcessorTargets
if (DEFINED PROC_FORCED_LABEL)
  set (PROC_LABEL ${PROC_FORCED_LABEL})
endif (DEFINED PROC_FORCED_LABEL)

# adding the proc flags to the build flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PROC_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PROC_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${PROC_FLAGS}")

if (WIN32)
    # we additionnally look in the MinGW path first then in the Gtkmm path,
    # so if you wish to build some dependant library, you have to install them in MinGW to use them
    set(CMAKE_PREFIX_PATH $ENV{MINGW_BASEPATH} $ENV{GTKMM_BASEPATH} CACHE STRING "Additional search paths")
endif (WIN32)

if (APPLE)
#    SET (CMAKE_OSX_ARCHITECTURES "i386;x86_64;" )
#    SET (CMAKE_TRY_COMPILE_OSX_ARCHITECTURES "i386;x86_64;" )
    SET (CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.5.sdk")
    SET (CMAKE_OSX_DEPLOYMENT_TARGET "10.5")
endif (APPLE)

option (AUTOMATED_BUILD_SYSTEM "TRUE if built by an automate" OFF)

option (BUILD_SHARED "Build rawtherapee with shared libraries" OFF)
option (WITH_BZIP "Build with Bzip2 support" ON)
option (WITH_MYFILE_MMAP "Build using memory mapped file" ON)
option (OPTION_OMP "Build with OpenMP support" ON)
option (BUILD_BUNDLE "Self-contained build" OFF)

# set install directories
if (WIN32 OR APPLE)
    set (BUILD_BUNDLE ON FORCE)
endif(WIN32 OR APPLE)

if (NOT DEFINED BUNDLE_BASE_INSTALL_DIR)
    set (BUNDLE_BASE_INSTALL_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/${CMAKE_BUILD_TYPE})
endif (NOT DEFINED BUNDLE_BASE_INSTALL_DIR)

if (BUILD_BUNDLE)
    set (BINDIR .)
    set (CMAKE_INSTALL_PREFIX ${BUNDLE_BASE_INSTALL_DIR})
endif (BUILD_BUNDLE)

if (NOT DEFINED BINDIR)
    set (BINDIR  ${CMAKE_INSTALL_PREFIX}/bin)
endif (NOT DEFINED BINDIR)

if (NOT DEFINED DATADIR)
    if (BUILD_BUNDLE)
        set (DATADIR .)
    else (BUILD_BUNDLE)
        set (DATADIR ${CMAKE_INSTALL_PREFIX}/share/rawtherapee)
    endif (BUILD_BUNDLE)
endif (NOT DEFINED DATADIR)

if (NOT DEFINED LIBDIR)
    if (BUILD_BUNDLE)
        set (LIBDIR .)
    else (BUILD_BUNDLE)
        # Respect CMAKE_INSTALL_LIBDIR if set
        if (DEFINED CMAKE_INSTALL_LIBDIR)
	    if (IS_ABSOLUTE ${LIBDIR})
            	set (LIBDIR ${CMAKE_INSTALL_LIBDIR})
	    else (IS_ABSOLUTE ${LIBDIR})
		set (LIBDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
	    endif (IS_ABSOLUTE ${LIBDIR})
        else (DEFINED CMAKE_INSTALL_LIBDIR)
            set (LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
        endif (DEFINED CMAKE_INSTALL_LIBDIR)
    endif (BUILD_BUNDLE)
endif (NOT DEFINED LIBDIR)

if (NOT DEFINED DOCDIR)
    if (BUILD_BUNDLE)
        set (DOCDIR ./doc)
    else (BUILD_BUNDLE)
        set (DOCDIR ${CMAKE_INSTALL_PREFIX}/share/doc/rawtherapee)
    endif (BUILD_BUNDLE)
endif (NOT DEFINED DOCDIR)

if (NOT DEFINED CREDITSDIR)
    if (BUILD_BUNDLE)
        set (CREDITSDIR .)
    else (BUILD_BUNDLE)
        set (CREDITSDIR ${CMAKE_INSTALL_PREFIX}/share/doc/rawtherapee)
    endif (BUILD_BUNDLE)
endif (NOT DEFINED CREDITSDIR)

if (NOT DEFINED LICENCEDIR)
    if (BUILD_BUNDLE)
        set (LICENCEDIR .)
    else (BUILD_BUNDLE)
        set (LICENCEDIR ${CMAKE_INSTALL_PREFIX}/share/doc/rawtherapee)
    endif (BUILD_BUNDLE)
endif (NOT DEFINED LICENCEDIR)

if (NOT DEFINED DESKTOPDIR)
    if (UNIX)
        set (DESKTOPDIR ${CMAKE_INSTALL_PREFIX}/share/applications)
    endif (UNIX)
endif (NOT DEFINED DESKTOPDIR)

if (NOT DEFINED ICONSDIR)
    if (UNIX)
        set (ICONSDIR ${CMAKE_INSTALL_PREFIX}/share/icons)
    endif (UNIX)
endif (NOT DEFINED ICONSDIR)

# non-bundle builds has to use absolute paths
if (NOT BUILD_BUNDLE AND
    NOT (IS_ABSOLUTE ${BINDIR} AND IS_ABSOLUTE ${DATADIR} AND IS_ABSOLUTE ${LIBDIR} AND
    IS_ABSOLUTE ${DOCDIR} AND IS_ABSOLUTE ${CREDITSDIR} AND IS_ABSOLUTE ${LICENCEDIR}))
    message (FATAL_ERROR "The paths has to be absolute or use -DBUILD_BUNDLE=ON")
endif ()

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
endif (WIN32)
# you may need lcms v1.xx for older version : pkg_check_modules (LCMS REQUIRED lcms<=1.99)
pkg_check_modules (LCMS REQUIRED lcms2)
pkg_check_modules (IPTCDATA REQUIRED libiptcdata)
find_package (JPEG REQUIRED)
find_package (PNG REQUIRED)
find_package (TIFF REQUIRED)
find_package (ZLIB REQUIRED)

if (AUTOMATED_BUILD_SYSTEM)

	message(WARNING "Automated build system mode! You have to make sure that AboutThisBuild.txt and rtgui/version.h are up-to-date by building RT at least one time with AUTOMATED_BUILD_SYSTEM=OFF.")

  # Checking for the presence of AboutThisBuild.txt
  find_file(ABOUT_FILE AboutThisBuild.txt PATHS "${CMAKE_CURRENT_SOURCE_DIR}" NO_DEFAULT_PATH)
  if (ABOUT_FILE STREQUAL ABOUT_FILE-NOTFOUND)
	  message(FATAL_ERROR "AboutThisBuild.txt not found! You have to build RT at least one time with AUTOMATED_BUILD_SYSTEM=OFF to generate this file.")
  endif (ABOUT_FILE STREQUAL ABOUT_FILE-NOTFOUND)

  # Checking for the presence of version.h
  find_file(VERSION_FILE version.h PATHS "${CMAKE_CURRENT_SOURCE_DIR}/rtgui" NO_DEFAULT_PATH)
  if (VERSION_FILE STREQUAL VERSION_FILE-NOTFOUND)
	  message(FATAL_ERROR "rtgui/version.h not found! You have to build RT at least one time with AUTOMATED_BUILD_SYSTEM=OFF to generate this file.")
  endif (VERSION_FILE STREQUAL VERSION_FILE-NOTFOUND)

else (AUTOMATED_BUILD_SYSTEM)

  # set the bit number information of the platform
  if (CMAKE_SIZEOF_VOID_P EQUAL 4)
    set(PROC_BIT_DEPTH 32 bits)
  elseif (CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(PROC_BIT_DEPTH 64 bits)
  endif (CMAKE_SIZEOF_VOID_P EQUAL 4)

  #generating AboutThisBuild.txt
  if (WIN32)

    add_custom_target (AboutFile ALL
      COMMAND ${CMAKE_COMMAND}
              -DSRC_DIR:STRING=${CMAKE_CURRENT_SOURCE_DIR}
              -DCACHE_NAME_SUFFIX:STRING=${CACHE_NAME_SUFFIX}
              -DSYSTEM:STRING=Windows
              -DPROC_LABEL:STRING=\"${PROC_LABEL}\"
              -DPROC_BIT_DEPTH:STRING=\"${PROC_BIT_DEPTH}\"
              -DBUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
              -DGTKMM_VERSION:STRING=${GTKMM_VERSION}
              -DCXX_FLAGS:STRING=\"${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${UPPER_CMAKE_BUILD_TYPE}}\"
              -DLFLAGS:STRING=\"${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_${UPPER_CMAKE_BUILD_TYPE}}\"
              -DOPTION_OMP:STRING=${OPTION_OMP}
              -DWITH_MYFILE_MMAP:STRING=${WITH_MYFILE_MMAP}
              -P ${CMAKE_CURRENT_SOURCE_DIR}/About-Windows.cmake
      COMMENT "Creating the about file"
    )
    find_file(HG_CMD hg.exe HINTS ENV Path PATH_SUFFIXES ../)

  elseif (APPLE)
 
    add_custom_target (AboutFile ALL
      COMMAND ${CMAKE_COMMAND}
              -DSRC_DIR:STRING=${CMAKE_CURRENT_SOURCE_DIR}
              -DCACHE_NAME_SUFFIX:STRING=${CACHE_NAME_SUFFIX}
              -DSYSTEM:STRING=Apple
              -DPROC_LABEL:STRING=\"${PROC_LABEL}\"
              -DPROC_BIT_DEPTH:STRING=\"${PROC_BIT_DEPTH}\"
              -DBUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
              -DGTKMM_VERSION:STRING=${GTKMM_VERSION}
              -DCXX_FLAGS:STRING=${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${UPPER_CMAKE_BUILD_TYPE}}
              -DLFLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_${UPPER_CMAKE_BUILD_TYPE}}
              -DOPTION_OMP:STRING=${OPTION_OMP}
              -DWITH_MYFILE_MMAP:STRING=${WITH_MYFILE_MMAP}
              -P ${CMAKE_CURRENT_SOURCE_DIR}/About-Apple.cmake
      COMMENT "Creating the about file"
    )
    find_file(HG_CMD hg PATHS "/opt/local/bin" "/usr/local/bin" "/usr/bin")
    find_file(HG_CMD hg)

  else (WIN32)
 
    add_custom_target (AboutFile ALL
      COMMAND ${CMAKE_COMMAND}
              -DSRC_DIR:STRING=${CMAKE_CURRENT_SOURCE_DIR}
              -DCACHE_NAME_SUFFIX:STRING=${CACHE_NAME_SUFFIX}
              -DSYSTEM:STRING=Linux
              -DPROC_LABEL:STRING=\"${PROC_LABEL}\"
              -DPROC_BIT_DEPTH:STRING=\"${PROC_BIT_DEPTH}\"
              -DBUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
              -DGTKMM_VERSION:STRING=${GTKMM_VERSION}
              -DCXX_FLAGS:STRING=${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${UPPER_CMAKE_BUILD_TYPE}}
              -DLFLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_${UPPER_CMAKE_BUILD_TYPE}}
              -DOPTION_OMP:STRING=${OPTION_OMP}
              -DWITH_MYFILE_MMAP:STRING=${WITH_MYFILE_MMAP}
              -P ${CMAKE_CURRENT_SOURCE_DIR}/About-Linux.cmake
      COMMENT "Creating the about file"
    )
    find_file(HG_CMD hg)

  endif (WIN32)
 
  add_dependencies(AboutFile Debug Release MinSizeRel RelWithDebInfo)

  if (HG_CMD STREQUAL HG_CMD-NOTFOUND)
    message(FATAL_ERROR "hg command not found!")
  else (HG_CMD STREQUAL HG_CMD-NOTFOUND)
    message(STATUS "hg command found: ${HG_CMD}")
  endif (HG_CMD STREQUAL HG_CMD-NOTFOUND)

endif (AUTOMATED_BUILD_SYSTEM)


# link witz bzip
if (WITH_BZIP)
  find_package(BZip2)
  if (BZIP2_FOUND)
    add_definitions (-DBZIP_SUPPORT)
    set (EXTRA_INCDIR ${EXTRA_LIB} ${BZIP2_INCLUDE_DIR})
    set (EXTRA_LIB ${EXTRA_LIB} ${BZIP2_LIBRARIES})
  endif (BZIP2_FOUND)
endif (WITH_BZIP)

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
install (FILES AboutThisBuild.txt DESTINATION ${CREDITSDIR})

add_subdirectory (rtexif)
add_subdirectory (rtengine)
add_subdirectory (rtgui)
add_subdirectory (rtdata)
