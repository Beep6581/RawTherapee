if(WIN32)
    cmake_minimum_required(VERSION 2.8.4)
else()
    cmake_minimum_required(VERSION 2.6)
endif()

# Must stay before the 'project' command:
if(${CMAKE_EXTRA_GENERATOR} MATCHES "Eclipse CDT4")
    set(CMAKE_CXX_COMPILER_ARG1 "-std=c++11" CACHE STRING "C++ version for eclipse" FORCE)
    # Users building with Eclipse should set CMAKE_ECLIPSE_VERSION through the
    # command line to their current version of Eclipse:
    #set(CMAKE_ECLIPSE_VERSION "4.6.0" CACHE STRING "Eclipse version" FORCE)
endif()

PROJECT(RawTherapee)

# The default target is Debug:
if(CMAKE_BUILD_TYPE STREQUAL "")
    set (CMAKE_BUILD_TYPE Debug CACHE STRING "One of: None Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()

string(TOUPPER ${CMAKE_BUILD_TYPE} UPPER_CMAKE_BUILD_TYPE)

# Set required C and C++ standards and check GCC version:
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.9")
    message(FATAL_ERROR "Building RawTherapee requires using GCC version 4.9 or higher!")
endif()

# We might want to build using the old C++ ABI, even when using a new GCC version:
if(USE_OLD_CXX_ABI)
    add_definitions (-D_GLIBCXX_USE_CXX11_ABI=0)
endif()

if(UPPER_CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    add_definitions (-D_DEBUG)
else()
    add_definitions (-DNDEBUG)
    add_definitions (-D_DNDEBUG)
endif()

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# Cache name suffix examples: "" = ~/.config/RawTherapee, "5" = ~/.config/RawTherapee-5, "_testing" = ~/.config/RawTherapee_testing
set(CACHE_NAME_SUFFIX "" CACHE STRING "RawTherapee's cache folder suffix")

# By default we don't use a specific processor target, so PROC_TARGET_NUMBER is set to 0.
# Specify other values to optimize for specific processor architecture as listed in ProcessorTargets.cmake:
set(PROC_TARGET_NUMBER 0 CACHE STRING "Selected target processor from the list above (taken from ProcessorTargets.cmake)")

# Set special compilation flags for rtengine which get added to CMAKE_CXX_FLAGS:
set(RTENGINE_CXX_FLAGS "" CACHE STRING "Special compilation flags for RTEngine")

# Loads the ProcessorTargets list:
include(ProcessorTargets.cmake)
set(PROC_LABEL "undefined" CACHE STRING "Target processor label, unused if PROC_TARGET_NUMBER = 0 or 2")
set(PROC_FLAGS "" CACHE STRING "Target processor related build/link flags")
if((NOT (PROC_TARGET_NUMBER EQUAL 0)) AND (NOT (PROC_TARGET_NUMBER EQUAL 2)))
    set(PROC_LABEL ${PROC_TARGET_${PROC_TARGET_NUMBER}_LABEL})
endif()
if(NOT(PROC_TARGET_NUMBER EQUAL 0))
    set(PROC_FLAGS ${PROC_TARGET_${PROC_TARGET_NUMBER}_FLAGS})
endif()
if(UNIX AND PROC_LABEL STREQUAL "undefined")
    execute_process(COMMAND uname -p OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE cpu)
    if("${cpu}" STREQUAL "unknown")
        set(PROC_LABEL "${CMAKE_SYSTEM_PROCESSOR}")
    else()
        set(PROC_LABEL "${cpu}")
    endif()
endif()

# If PROC_FORCED_LABEL exists, its value is loaded in PROC_LABEL to override the one from ProcessorTargets:
if(DEFINED PROC_FORCED_LABEL)
    set(PROC_LABEL ${PROC_FORCED_LABEL})
endif()

# Add the proc flags to the build flags:
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PROC_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PROC_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${PROC_FLAGS}")

# Stop compilation on typos such as std:swap (missing colon will be detected as unused label):
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=unused-label")

if(WIN32)
    # Add additional paths. Look in the MinGW path first, then in the Gtkmm path.
    # If you wish to build some dependent libraries, you have to install them in MinGW to use them:
    set(CMAKE_PREFIX_PATH $ENV{MINGW_BASEPATH} $ENV{GTKMM_BASEPATH} CACHE STRING "Additional search paths")
endif()

if(APPLE)
    if(CMAKE_CXX_COMPILER MATCHES "g\\+\\+-mp-4.[5-8]" OR CMAKE_CXX_COMPILER_ARG1 MATCHES "g\\+\\+-mp-4.[5-8]")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /usr/lib/libstdc++.6.dylib")
        message(STATUS "CMAKE_CXX_COMPILER is MacPorts GCC.\n   CMAKE_EXE_LINKER_FLAGS: ${CMAKE_EXE_LINKER_FLAGS}")
    endif()
endif()

option(USE_EXPERIMENTAL_LANG_VERSIONS "Build RT with -std=c++0x" OFF)
option(BUILD_SHARED "Build rawtherapee with shared libraries" OFF)
option(WITH_BZIP "Build with Bzip2 support" ON)
option(WITH_MYFILE_MMAP "Build using memory mapped file" ON)
option(WITH_LTO "Build with link-time optimizations" OFF)
option(WITH_SAN "Build with run-time sanitizer" OFF)
option(WITH_PROF "Build with profiling instrumentation" OFF)
option(OPTION_OMP "Build with OpenMP support" ON)
option(STRICT_MUTEX "True (recommended): MyMutex will behave like POSIX Mutex; False: MyMutex will behave like POSIX RecMutex; Note: forced to ON for Debug builds" ON)
option(TRACE_MYRWMUTEX "Trace RT's custom R/W Mutex (Debug builds only); redirecting std::out to a file is strongly recommended!" OFF)
option(AUTO_GDK_FLUSH "Use gdk_flush on all gdk_thread_leave other than the GUI thread; set it ON if you experience X Server warning/errors" OFF)

# Set installation directories:
if(WIN32 OR APPLE)
    if(BUILD_BUNDLE)
        message(STATUS "You have set BUILD_BUNDLE=ON but this is not necessary - the option if force to ON on Windows and macOS.")
    endif()
    set (BUILD_BUNDLE ON FORCE)
endif()

if(NOT DEFINED BUNDLE_BASE_INSTALL_DIR)
    if(APPLE)
        set(BUNDLE_BASE_INSTALL_DIR "${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/MacOS")
    else()
        set(BUNDLE_BASE_INSTALL_DIR "${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}")
    endif()
endif()

if(BUILD_BUNDLE)
    set(BINDIR .)
    set(CMAKE_INSTALL_PREFIX "${BUNDLE_BASE_INSTALL_DIR}")
endif()

if(NOT DEFINED BINDIR)
    set(BINDIR "${CMAKE_INSTALL_PREFIX}/bin")
endif()

if(NOT DEFINED DATADIR)
    if(BUILD_BUNDLE)
        if(APPLE)
            set(DATADIR "../Resources")
        else()
            set(DATADIR .)
        endif()
    else()
        set(DATADIR "${CMAKE_INSTALL_PREFIX}/share/rawtherapee")
    endif()
endif()

if(NOT DEFINED LIBDIR)
    if(BUILD_BUNDLE)
        if(APPLE)
            set(LIBDIR "../Frameworks")
        else()
            set(LIBDIR .)
        endif()
    else()
        # Respect CMAKE_INSTALL_LIBDIR if set
        if(DEFINED CMAKE_INSTALL_LIBDIR)
            if(IS_ABSOLUTE "${LIBDIR}")
                set(LIBDIR "${CMAKE_INSTALL_LIBDIR}")
            else()
                set(LIBDIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
            endif()
        else()
            set(LIBDIR "${CMAKE_INSTALL_PREFIX}/lib")
        endif()
    endif()
endif()

if(NOT DEFINED DOCDIR)
    if(BUILD_BUNDLE)
        set(DOCDIR "${DATADIR}/share/doc")
    else()
        set(DOCDIR "${CMAKE_INSTALL_PREFIX}/share/doc/rawtherapee")
    endif()
endif()

if(NOT DEFINED CREDITSDIR)
    if(BUILD_BUNDLE)
        set(CREDITSDIR "${DATADIR}")
    else()
        set(CREDITSDIR "${CMAKE_INSTALL_PREFIX}/share/doc/rawtherapee")
    endif()
endif()

if(NOT DEFINED LICENCEDIR)
    if(BUILD_BUNDLE)
        set(LICENCEDIR "${DATADIR}")
    else()
        set(LICENCEDIR "${CMAKE_INSTALL_PREFIX}/share/doc/rawtherapee")
    endif()
endif()

if(NOT DEFINED DESKTOPDIR)
    if(UNIX)
        set(DESKTOPDIR "${CMAKE_INSTALL_PREFIX}/share/applications") # Absolute
        #set(DESKTOPDIR "${DATADIR}/share/applications") # Relative
    endif()
endif()

if(NOT DEFINED ICONSDIR)
    if(UNIX)
        set(ICONSDIR "${CMAKE_INSTALL_PREFIX}/share/icons")
    endif()
endif()

if(NOT DEFINED APPDATADIR)
    if(UNIX)
        set(APPDATADIR "${CMAKE_INSTALL_PREFIX}/share/appdata")
    endif()
endif()

# Enforce absolute paths for non-bundle builds:
if(NOT BUILD_BUNDLE)
    if(NOT(IS_ABSOLUTE "${BINDIR}"))
        message(FATAL_ERROR "The BINDIR path has to be absolute, or use -DBUILD_BUNDLE=ON")
    elseif(NOT(IS_ABSOLUTE "${DATADIR}"))
        message(FATAL_ERROR "The DATADIR path has to be absolute, or use -DBUILD_BUNDLE=ON")
    elseif(NOT(IS_ABSOLUTE "${LIBDIR}"))
        message(FATAL_ERROR "The LIBDIR path has to be absolute, or use -DBUILD_BUNDLE=ON")
    elseif(NOT(IS_ABSOLUTE "${DOCDIR}"))
        message(FATAL_ERROR "The DOCDIR path has to be absolute, or use -DBUILD_BUNDLE=ON")
    elseif(NOT(IS_ABSOLUTE "${CREDITSDIR}"))
        message(FATAL_ERROR "The CREDITSDIR path has to be absolute, or use -DBUILD_BUNDLE=ON")
    elseif(NOT(IS_ABSOLUTE "${LICENCEDIR}"))
        message(FATAL_ERROR "The LICENCEDIR path has to be absolute, or use -DBUILD_BUNDLE=ON")
    endif()
endif()

# MyMutex:
if(STRICT_MUTEX OR UPPER_CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    add_definitions(-DSTRICT_MUTEX=1)
else()
    add_definitions(-DSTRICT_MUTEX=0)
endif()

# MyRWMutex:
if(TRACE_MYRWMUTEX)
    add_definitions(-DTRACE_MYRWMUTEX=1)
else()
    add_definitions(-DTRACE_MYRWMUTEX=0)
endif()

if(AUTO_GDK_FLUSH)
    add_definitions(-DAUTO_GDK_FLUSH=1)
else()
    add_definitions(-DAUTO_GDK_FLUSH=0)
endif()

# Check for libraries:
find_package(PkgConfig)
pkg_check_modules (GTK     REQUIRED gtk+-3.0>=3.16)
pkg_check_modules (GLIB2   REQUIRED glib-2.0>=2.44)
pkg_check_modules (GLIBMM  REQUIRED glibmm-2.4>=2.44)
pkg_check_modules (GTKMM   REQUIRED gtkmm-3.0>=3.16)
pkg_check_modules (CAIROMM REQUIRED cairomm-1.0)
pkg_check_modules (GIO     REQUIRED gio-2.0>=2.44)
pkg_check_modules (GIOMM   REQUIRED giomm-2.4>=2.44)
pkg_check_modules (GTHREAD REQUIRED gthread-2.0>=2.44)
pkg_check_modules (GOBJECT REQUIRED gobject-2.0>=2.44)
pkg_check_modules (SIGC    REQUIRED sigc++-2.0>=2.3.1)

if(WIN32)
    add_definitions(-DWIN32)
    add_definitions(-D_WIN32)
    if(MINGW)
        add_definitions(-D__MINGW32__)
    endif()
    if(CMAKE_SIZEOF_VOID_P EQUAL 4)
        add_definitions(-DWINVER=0x0501)
    endif()
    set(EXTRA_LIB "-lws2_32 -lshlwapi")
endif()

pkg_check_modules(LCMS REQUIRED lcms2>=2.6)
pkg_check_modules(EXPAT REQUIRED expat>=2.1)
pkg_check_modules(FFTW3F REQUIRED fftw3f)
pkg_check_modules(IPTCDATA REQUIRED libiptcdata)
find_package(JPEG REQUIRED)
find_package(PNG REQUIRED)
find_package(TIFF REQUIRED)
find_package(ZLIB REQUIRED)

# Link with bzip:
if(WITH_BZIP)
    find_package(BZip2)
    if(BZIP2_FOUND)
        add_definitions(-DBZIP_SUPPORT)
        set(EXTRA_INCDIR ${BZIP2_INCLUDE_DIR})
        set(EXTRA_LIB ${EXTRA_LIB} ${BZIP2_LIBRARIES})
    endif()
endif()

# Check for libcanberra-gtk3 (sound events on Linux):
if(UNIX AND(NOT APPLE))
    pkg_check_modules(CANBERRA-GTK REQUIRED libcanberra-gtk3)
endif()

if(WITH_MYFILE_MMAP)
    add_definitions(-DMYFILE_MMAP)
endif()

if(WITH_LTO)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -flto")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto")
endif()

if(WITH_SAN)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=${WITH_SAN}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=${WITH_SAN}")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=${WITH_SAN}")
endif()

if(WITH_PROF)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
endif()

if(OPTION_OMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Werror=unknown-pragmas -Wall -Wno-unused-result -Wno-deprecated-declarations")
    endif()
endif()

# Find out whether we are building out of source:
get_filename_component(ABS_SOURCE_DIR "${PROJECT_SOURCE_DIR}" ABSOLUTE)
get_filename_component(ABS_BINARY_DIR "${CMAKE_BINARY_DIR}" ABSOLUTE)
set(OUT_OF_SOURCE_BUILD TRUE)
if(ABS_SOURCE_DIR STREQUAL ABS_BINARY_DIR)
    set(OUT_OF_SOURCE_BUILD FALSE)
    message(WARNING "You are performing an in-source build. This is discouraged. For an explanation and the advantages of out-of-source builds, please refer to http://www.cmake.org/Wiki/CMake_FAQ#What_is_an_.22out-of-source.22_build.3F")
endif()

# Remove files which could require manual work.
# We will remove this after some time to have a clean build system without file modifications in the source tree again. (?)
set(OOSB_FILES "${PROJECT_SOURCE_DIR}/rtdata/rawtherapee.desktop" "${PROJECT_SOURCE_DIR}/rtgui/version.h" "${PROJECT_SOURCE_DIR}/rtgui/config.h" "${PROJECT_SOURCE_DIR}/AboutThisBuild.txt")
if(OUT_OF_SOURCE_BUILD)
    foreach(f ${OOSB_FILES})
        file(REMOVE "${f}")
    endforeach()
endif()

# Check for generated files in the source tree which should not be there when doing an out-of-source build.
# Without checking for this it might happen that old versions are used for the compilation:
if(OUT_OF_SOURCE_BUILD)
    foreach(f ${OOSB_FILES})
        if(EXISTS "${f}")
            message(SEND_ERROR "Generated \"${f}\" found inside the source tree. Please remove it as it is a relic of the old build system and prevents valid compilation now.")
        endif()
    endforeach()
endif()

### Start generating AboutThisBuild.txt
# Set the platform bit-depth:
if(CMAKE_SIZEOF_VOID_P EQUAL 4)
    set(PROC_BIT_DEPTH 32 bits)
elseif()
    set(PROC_BIT_DEPTH 64 bits)
endif()

# Get compiler name and version.
# Only CMake > 2.8.7 knows CMAKE_*_COMPILER_VERSION
if(CMAKE_VERSION VERSION_GREATER 2.8.7)
    get_filename_component(COMPILER_INFO ${CMAKE_C_COMPILER} NAME_WE)
    set(COMPILER_INFO "${COMPILER_INFO} ${CMAKE_C_COMPILER_VERSION}")
else()
    execute_process(COMMAND gcc -dumpversion OUTPUT_VARIABLE GCC_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(COMPILER_INFO "gcc ${GCC_VERSION}")
endif()

# Get C++ and linker flags for rtengine (the GUI's C++ flags may have fewer flags):
set(CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${UPPER_CMAKE_BUILD_TYPE}} ${RTENGINE_CXX_FLAGS}")
set(LFLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_${UPPER_CMAKE_BUILD_TYPE}}")

set(ABOUT_COMMAND_WITH_ARGS ${CMAKE_COMMAND}
    -DPROJECT_SOURCE_DIR:STRING=${PROJECT_SOURCE_DIR}
    -DCACHE_NAME_SUFFIX:STRING=${CACHE_NAME_SUFFIX}
    -DPROC_LABEL:STRING="${PROC_LABEL}"
    -DPROC_BIT_DEPTH:STRING="${PROC_BIT_DEPTH}"
    -DBUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DGTKMM_VERSION:STRING=${GTKMM_VERSION}
    -DOPTION_OMP:STRING=${OPTION_OMP}
    -DWITH_MYFILE_MMAP:STRING=${WITH_MYFILE_MMAP})

if(WIN32)
    list(APPEND ABOUT_COMMAND_WITH_ARGS -DSYSTEM:STRING=Windows
        -DCXX_FLAGS:STRING="${CXX_FLAGS}"
        -DLFLAGS:STRING="${LFLAGS}"
        -DCOMPILER_INFO:STRING="${COMPILER_INFO}"
        -DCMAKE_INSTALL_PREFIX:STRING="${CMAKE_INSTALL_PREFIX}"
        -DBIT_DEPTH:STRING="${CMAKE_SIZEOF_VOID_P}")
elseif(APPLE)
    list(APPEND ABOUT_COMMAND_WITH_ARGS -DSYSTEM:STRING=Apple
        -DCXX_FLAGS:STRING=${CXX_FLAGS}
        -DLFLAGS:STRING=${LFLAGS}
        -DCOMPILER_INFO:STRING=${COMPILER_INFO})
else()
    list(APPEND ABOUT_COMMAND_WITH_ARGS -DSYSTEM:STRING=Linux
        -DCXX_FLAGS:STRING=${CXX_FLAGS}
        -DLFLAGS:STRING=${LFLAGS}
        -DCOMPILER_INFO:STRING=${COMPILER_INFO})
endif()

list(APPEND ABOUT_COMMAND_WITH_ARGS -P "${PROJECT_SOURCE_DIR}/UpdateInfo.cmake")

add_custom_target(UpdateInfo ALL
    COMMAND ${ABOUT_COMMAND_WITH_ARGS}
    COMMENT "Creating AboutThisBuild.txt and other version-dependent files")
### End generating AboutThisBuild.txt

install(FILES AUTHORS.txt DESTINATION "${CREDITSDIR}")
install(FILES LICENSE.txt DESTINATION "${LICENCEDIR}")
install(FILES "${CMAKE_BINARY_DIR}/AboutThisBuild.txt" DESTINATION "${CREDITSDIR}")
install(FILES RELEASE_NOTES.txt DESTINATION "${CREDITSDIR}" OPTIONAL)

if(UNIX OR WIN32)
    install(FILES "${PROJECT_SOURCE_DIR}/doc/manpage/rawtherapee.1" DESTINATION "${DOCDIR}/man")
endif()

if(WIN32)
    install(DIRECTORY "licenses" DESTINATION "${LICENCEDIR}")
endif()

if(UNIX)
    install(FILES rawtherapee.appdata.xml DESTINATION "${APPDATADIR}")
endif()

add_subdirectory(rtexif)
add_subdirectory(rtengine)
add_subdirectory(rtgui)
add_subdirectory(rtdata)

