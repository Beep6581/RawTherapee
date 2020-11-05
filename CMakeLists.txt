# Use of SOURCE_DIR target property begins at cmake 3.7
if(APPLE)
    cmake_minimum_required(VERSION 3.7)
else()
    cmake_minimum_required(VERSION 3.5)
endif()

# Compiler id for Apple Clang
if(APPLE)
    cmake_policy(SET CMP0025 NEW)
    cmake_policy(SET CMP0037 NEW)
endif()

# Must stay before the project() command:
if(${CMAKE_EXTRA_GENERATOR} MATCHES "Eclipse CDT4")
    set(CMAKE_CXX_COMPILER_ARG1
        "-std=c++11"
        CACHE STRING "C++ version for eclipse" FORCE)
    # Users building with Eclipse should set CMAKE_ECLIPSE_VERSION through the
    # command line to their current version of Eclipse:
    # set(CMAKE_ECLIPSE_VERSION "4.6.0" CACHE STRING "Eclipse version" FORCE)
endif()

project(RawTherapee)

# The default target is Debug:
if(CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE
        Debug
        CACHE STRING "One of: None Debug Release RelWithDebInfo MinSizeRel"
              FORCE)
endif()

string(TOUPPER ${CMAKE_BUILD_TYPE} UPPER_CMAKE_BUILD_TYPE)

# Set required C and C++ standards and check GCC version:
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION
                                            VERSION_LESS "4.9")
    message(
        FATAL_ERROR
            "Building RawTherapee requires using GCC version 4.9 or higher!")
endif()

# Warning for GCC 10, which causes problems #5749:
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL "10.1")
    message(STATUS "WARNING: gcc ${CMAKE_CXX_COMPILER_VERSION} is known to miscompile RawTherapee when using -ftree-loop-vectorize, forcing the option to be off")    
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-tree-loop-vectorize")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-tree-loop-vectorize")
endif()

# We might want to build using the old C++ ABI, even when using a new GCC
# version:
if(USE_OLD_CXX_ABI)
    add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)
endif()

if(UPPER_CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    add_definitions(-D_DEBUG)
else()
    add_definitions(-DNDEBUG)
    add_definitions(-D_DNDEBUG)
endif()

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# Cache name suffix examples: "" = ~/.config/RawTherapee, "5" =
# ~/.config/RawTherapee-5, "_testing" = ~/.config/RawTherapee_testing Use "" for
# stable releases and "5-dev" for anything else.
set(CACHE_NAME_SUFFIX
    ""
    CACHE STRING "RawTherapee's cache folder suffix")

# For macOS only, OSX_DEV_BUILD option allows using relative paths instead of absolute
# paths. Consequently, for development builds, application can be launching without
# being bundled. However, file access can be restricted for some folder.
option(OSX_DEV_BUILD "Generate macOS development builds" OFF)

# By default we don't use a specific processor target, so PROC_TARGET_NUMBER is
# set to 0. Specify other values to optimize for specific processor architecture
# as listed in ProcessorTargets.cmake:
set(PROC_TARGET_NUMBER
    0
    CACHE
        STRING
        "Selected target processor from the list above (taken from ProcessorTargets.cmake)"
)

# Set special compilation flags for rtengine which get added to CMAKE_CXX_FLAGS:
# Some Linux distros build with -O2 instead of -O3. We explicitly enable auto
# vectorization by using -ftree-vectorize
set(RTENGINE_CXX_FLAGS
    "-ftree-vectorize"
    CACHE STRING "Special compilation flags for RTEngine")

# Loads the ProcessorTargets list:
include(ProcessorTargets.cmake)
set(PROC_LABEL
    "undefined"
    CACHE STRING
          "Target processor label, unused if PROC_TARGET_NUMBER = 0 or 2")
set(PROC_FLAGS
    ""
    CACHE STRING "Target processor related build/link flags")
if((NOT (PROC_TARGET_NUMBER EQUAL 0)) AND (NOT (PROC_TARGET_NUMBER EQUAL 2)))
    set(PROC_LABEL ${PROC_TARGET_${PROC_TARGET_NUMBER}_LABEL})
endif()
if(NOT (PROC_TARGET_NUMBER EQUAL 0))
    set(PROC_FLAGS ${PROC_TARGET_${PROC_TARGET_NUMBER}_FLAGS})
endif()
if(UNIX AND PROC_LABEL STREQUAL "undefined")
    execute_process(
        COMMAND uname -p
        OUTPUT_STRIP_TRAILING_WHITESPACE
        OUTPUT_VARIABLE cpu)
    if("${cpu}" STREQUAL "unknown")
        set(PROC_LABEL "${CMAKE_SYSTEM_PROCESSOR}")
    else()
        set(PROC_LABEL "${cpu}")
    endif()
endif()

# If PROC_FORCED_LABEL exists, its value is loaded in PROC_LABEL to override the
# one from ProcessorTargets:
if(DEFINED PROC_FORCED_LABEL)
    set(PROC_LABEL ${PROC_FORCED_LABEL})
endif()

# Add the proc flags to the build flags:
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PROC_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PROC_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${PROC_FLAGS}")

# Stop compilation on typos such as std:swap (missing colon will be detected as
# unused label):
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=unused-label")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=delete-incomplete")

# Do net set math errno, as we never check its value.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-math-errno")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-math-errno")

# Special treatment for x87 and x86-32 SSE (see GitHub issue #4324)
include(FindX87Math)
if(HAVE_X87_MATH)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffloat-store")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffloat-store")
endif()
if(HAVE_X86_SSE_MATH)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2 -mfpmath=sse")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse2 -mfpmath=sse")
endif()

# On i386 Linux we can fix unaligned SSE malloc (see GitHub issue #4432)
include(FindUnalignedMalloc)

if(WIN32)
    # Add additional paths. Look in the MinGW path first, then in the Gtkmm
    # path. If you wish to build some dependent libraries, you have to install
    # them in MinGW to use them:
    set(CMAKE_PREFIX_PATH
        $ENV{MINGW_BASEPATH} $ENV{GTKMM_BASEPATH}
        CACHE STRING "Additional search paths")
endif()

if(APPLE) # TODO make -mtune generic conditional and/or specifiable.
    if(CMAKE_CXX_COMPILER MATCHES "g\\+\\+-mp-4.[5-8]"
       OR CMAKE_CXX_COMPILER_ARG1 MATCHES "g\\+\\+-mp-4.[5-8]")
        set(CMAKE_EXE_LINKER_FLAGS
            "${CMAKE_EXE_LINKER_FLAGS} /usr/lib/libstdc++.6.dylib -Wl,-headerpad_max_install_names -mtune=generic"
        )
        message(
            STATUS
                "CMAKE_CXX_COMPILER is MacPorts GCC.\n   CMAKE_EXE_LINKER_FLAGS: ${CMAKE_EXE_LINKER_FLAGS}"
        )
    endif()

    # Set minimum system version
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -mmacosx-version-min=${CMAKE_OSX_DEPLOYMENT_TARGET}"
    )
    set(CMAKE_C_FLAGS
        "${CMAKE_C_FLAGS} -mmacosx-version-min=${CMAKE_OSX_DEPLOYMENT_TARGET}")
    set(CMAKE_EXE_LINKER_FLAGS
        "${CMAKE_EXE_LINKER_FLAGS} -Wl,-headerpad_max_install_names -mtune=generic"
    )
endif()

option(USE_EXPERIMENTAL_LANG_VERSIONS "Build with -std=c++0x" OFF)
option(BUILD_SHARED "Build with shared libraries" OFF)
option(WITH_BENCHMARK "Build with benchmark code" OFF)
option(WITH_MYFILE_MMAP "Build using memory mapped file" ON)
option(WITH_LTO "Build with link-time optimizations" OFF)
option(WITH_SAN "Build with run-time sanitizer" OFF)
option(WITH_PROF "Build with profiling instrumentation" OFF)
option(WITH_SYSTEM_KLT "Build using system KLT library." OFF)
option(OPTION_OMP "Build with OpenMP support" ON)
option(
    STRICT_MUTEX
    "True (recommended): MyMutex will behave like POSIX Mutex; False: MyMutex will behave like POSIX RecMutex; Note: forced to ON for Debug builds"
    ON)
option(
    TRACE_MYRWMUTEX
    "Trace custom R/W Mutex (Debug builds only); redirecting std::out to a file is strongly recommended!"
    OFF)
option(
    AUTO_GDK_FLUSH
    "Use gdk_flush on all gdk_thread_leave other than the GUI thread; set it ON if you experience X Server warning/errors"
    OFF)
# option(TARGET32BIT "Build for 32-bit architecture when ON, otherwise 64-bit.
# Default is OFF" OFF)

option(ENABLE_TCMALLOC "Use the tcmalloc library if available" OFF)
set(TCMALLOC_LIB_DIR
    ""
    CACHE PATH "Custom path for the tcmalloc library")

# Set installation directories:
if(WIN32 OR APPLE)
    if(BUILD_BUNDLE)
        message(
            STATUS
                "You have set BUILD_BUNDLE=ON but this is not necessary - the option is forced to ON for Windows and macOS."
        )
    endif()
    set(BUILD_BUNDLE ON FORCE)
endif()

if(NOT DEFINED BUNDLE_BASE_INSTALL_DIR)
    if(APPLE)
        set(BUNDLE_BASE_INSTALL_DIR "${PROJECT_SOURCE_DIR}/build/${CMAKE_BUILD_TYPE}/MacOS")
    else()
        set(BUNDLE_BASE_INSTALL_DIR "${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}")
    endif()
endif()

if(BUILD_BUNDLE)
    if(APPLE)
        set(BINDIR "${BUNDLE_BASE_INSTALL_DIR}")
        set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/build")
    else()
        set(BINDIR .)
        set(CMAKE_INSTALL_PREFIX "${BUNDLE_BASE_INSTALL_DIR}")
    endif()
endif()

if(NOT DEFINED BINDIR)
    set(BINDIR "${CMAKE_INSTALL_PREFIX}/bin")
endif()

if(NOT DEFINED DATADIR)
    if(BUILD_BUNDLE)
        if(APPLE)
            set(DATADIR "${PROJECT_SOURCE_DIR}/build/${CMAKE_BUILD_TYPE}/Resources/share")
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
            set(LIBDIR "${PROJECT_SOURCE_DIR}/build/${CMAKE_BUILD_TYPE}/Frameworks")
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
        if(APPLE)
            set(DOCDIR "${DATADIR}/doc")
        else()
            set(DOCDIR "${DATADIR}/share/doc")
        endif()
    else()
        set(DOCDIR "${CMAKE_INSTALL_PREFIX}/share/doc/rawtherapee")
    endif()
endif()

if(NOT DEFINED CREDITSDIR)
    if(BUILD_BUNDLE)
        if(APPLE)
            set(CREDITSDIR "${DATADIR}")
        else()
            set(CREDITSDIR "${DATADIR}")
        endif()
    else()
        set(CREDITSDIR "${CMAKE_INSTALL_PREFIX}/share/doc/rawtherapee")
    endif()
endif()

if(NOT DEFINED LICENCEDIR)
    if(BUILD_BUNDLE)
        if(APPLE)
            set(LICENCEDIR "${CREDITSDIR}")
        else()
            set(LICENCEDIR "${DATADIR}")
        endif()
    else()
        set(LICENCEDIR "${CMAKE_INSTALL_PREFIX}/share/doc/rawtherapee")
    endif()
endif()

if(NOT DEFINED DESKTOPDIR)
    if(UNIX)
        if(BUILD_BUNDLE)
            if(APPLE)
                set(DESKTOPDIR "${DATADIR}/applications")
            else()
                set(DESKTOPDIR "${DATADIR}/share/applications")
            endif()
        else()
            set(DESKTOPDIR "${CMAKE_INSTALL_PREFIX}/share/applications")
        endif()
    endif()
endif()

if(NOT DEFINED ICONSDIR)
    if(UNIX)
        if(BUILD_BUNDLE)
            if(APPLE)
                set(ICONSDIR "${DATADIR}/icons")
            else()
                set(ICONSDIR "${DATADIR}/share/icons")
            endif()
        else()
            set(ICONSDIR "${CMAKE_INSTALL_PREFIX}/share/icons")
        endif()
    endif()
endif()

if(NOT DEFINED APPDATADIR)
    if(UNIX)
        if(BUILD_BUNDLE)
            if(APPLE)
                set(APPDATADIR "${DATADIR}/metainfo")
            else()
                set(APPDATADIR "${DATADIR}/share/metainfo")
            endif()
        else()
            set(APPDATADIR "${CMAKE_INSTALL_PREFIX}/share/metainfo")
        endif()
    endif()
endif()

if(NOT APPLE)
    if(DEFINED LENSFUNDBDIR AND NOT IS_ABSOLUTE "${LENSFUNDBDIR}")
        set(LENSFUNDBDIR "${DATADIR}/${LENSFUNDBDIR}")
    endif()
endif()

if(APPLE)
    if("${CODESIGNID}")
        set(CODESIGNID
            "${CODESIGNID}"
            CACHE STRING "Codesigning Identity")
    endif()
    if("${NOTARY}")
        set(NOTARY
            "${NOTARY}"
            CACHE STRING "Notarization Identity")
    endif()
    if("${LOCAL_PREFIX}")
        set(LOCAL_PREFIX
            "${LOCAL_PREFIX}"
            CACHE STRING "macos/gtk parent directory ie /usr or /opt")
    else()
        set(LOCAL_PREFIX
            /usr
            CACHE STRING "macos/gtk parent directory ie /usr or /opt")
    endif()
    if("${FANCY_DMG}")
        set(FANCY_DMG
            ON
            CACHE BOOL
                  "Use the andreyvit/create-dmg script to make a fancy .dmg")
    endif()
endif()

# Enforce absolute paths for non-bundle builds:
if(NOT BUILD_BUNDLE)
    foreach(path BINDIR DATADIR LIBDIR DOCDIR CREDITSDIR LICENCEDIR)
        if(NOT (IS_ABSOLUTE "${${path}}"))
            message(
                FATAL_ERROR
                    "The ${path} path has to be absolute when using -DBUILD_BUNDLE=OFF"
            )
        endif()
    endforeach()
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

if(WIN32)
    pkg_check_modules(GTK REQUIRED gtk+-3.0>=3.22.24)
    pkg_check_modules(GTKMM REQUIRED gtkmm-3.0>=3.22)
else()
    pkg_check_modules(GTK REQUIRED gtk+-3.0>=3.16)
    pkg_check_modules(GTKMM REQUIRED gtkmm-3.0>=3.16)
endif()

if(GTK_VERSION VERSION_GREATER "3.24.1" AND GTK_VERSION VERSION_LESS "3.24.7")
    if(GTK_VERSION VERSION_EQUAL "3.24.5")
        set(CERTAINTY "known to")
    else()
        set(CERTAINTY "likely to")
    endif()
    message(
        WARNING
            "\nWarning! You are using GTK+ version "
            ${GTK_VERSION}
            " which is "
            ${CERTAINTY}
            " have an issue where combobox menu scroll-arrows are missing when a Gtk::ComboBox list does not fit vertically on the screen. As a result, users of your build will not be able to select items in the following comboboxes: Processing Profiles, Film Simulation, and the camera and lens profiles in Profiled Lens Correction.\nIt is recommended that you either downgrade GTK+ to <= 3.24.1 or upgrade to >= 3.24.7."
    )
endif()

pkg_check_modules(GLIB2 REQUIRED glib-2.0>=2.44)
pkg_check_modules(GLIBMM REQUIRED glibmm-2.4>=2.44)
pkg_check_modules(CAIROMM REQUIRED cairomm-1.0)
pkg_check_modules(GIO REQUIRED gio-2.0>=2.44)
pkg_check_modules(GIOMM REQUIRED giomm-2.4>=2.44)
pkg_check_modules(GTHREAD REQUIRED gthread-2.0>=2.44)
pkg_check_modules(GOBJECT REQUIRED gobject-2.0>=2.44)
pkg_check_modules(SIGC REQUIRED sigc++-2.0>=2.3.1)
pkg_check_modules(LENSFUN REQUIRED lensfun>=0.2)
pkg_check_modules(RSVG REQUIRED librsvg-2.0>=2.40)

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
pkg_check_modules(TIFF REQUIRED libtiff-4>=4.0.4)
find_package(JPEG REQUIRED)
find_package(PNG REQUIRED)
find_package(ZLIB REQUIRED)
if(WITH_SYSTEM_KLT)
    find_package(KLT REQUIRED)
endif()

# Check for libcanberra-gtk3 (sound events on Linux):
if(UNIX AND (NOT APPLE))
    option(USE_LIBCANBERRA "Build with libcanberra" ON)
    if(USE_LIBCANBERRA)
       pkg_check_modules(CANBERRA-GTK REQUIRED libcanberra-gtk3)
       add_definitions(-DUSE_CANBERRA)
    endif()
endif()

if(WITH_MYFILE_MMAP)
    add_definitions(-DMYFILE_MMAP)
endif()

if(WITH_LTO)
    # Using LTO with older versions of binutils requires setting extra flags
    set(BINUTILS_VERSION_MININUM "2.29")
    execute_process(COMMAND ar --version
                    OUTPUT_VARIABLE BINUTILS_VERSION_DETECTED)
    string(REGEX REPLACE ".* ([0-9.]+)\n.*" "\\1" BINUTILS_VERSION_DETECTED
                         "${BINUTILS_VERSION_DETECTED}")
    if("${BINUTILS_VERSION_DETECTED}" VERSION_LESS
       "${BINUTILS_VERSION_MININUM}")
        if(APPLE)
            if(!CMAKE_AR)
                set(CMAKE_AR "/opt/local/bin/ar")
            endif()
            if(!CMAKE_RANLIB)
                set(CMAKE_RANLIB "/opt/local/bin/ranlib")
            endif()
        else()
            set(CMAKE_AR "/usr/bin/gcc-ar")
            set(CMAKE_RANLIB "/usr/bin/gcc-ranlib")
        endif()
        message(
            STATUS
                "Binutils version detected as less than "
                ${BINUTILS_VERSION_MININUM}
                " - setting CMake parameters to enable LTO linking:\n    CMAKE_AR=\""
                ${CMAKE_AR}
                "\"\n    CMAKE_RANLIB=\""
                ${CMAKE_RANLIB}
                "\"")
    endif()

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -flto")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto")
endif()

if(WITH_SAN)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=${WITH_SAN}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=${WITH_SAN}")
    set(CMAKE_EXE_LINKER_FLAGS
        "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=${WITH_SAN}")
endif()

if(WITH_PROF)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
endif()

set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Wall -Wuninitialized -Wcast-qual -Wno-deprecated-declarations -Wno-unused-result -Wunused-macros"
)
if(OPTION_OMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Werror=unknown-pragmas")
    endif()
endif()

# check for libfftw3f_omp
include(CheckCSourceCompiles)
if(OPENMP_FOUND)
    set(CMAKE_REQUIRED_INCLUDES ${FFTW3F_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LIBRARIES)
    foreach(l ${FFTW3F_LIBRARIES})
        find_library(_f ${l} PATHS ${FFTW3F_LIBRARY_DIRS})
        set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${_f})
    endforeach()
    check_c_source_compiles(
        "#include <fftw3.h>
int main()
{
    fftwf_init_threads();
    fftwf_plan_with_nthreads(1);
    return 0;
}"
        _fftw3f_multithread)
    if(_fftw3f_multithread)
        add_definitions(-DRT_FFTW3F_OMP)
    else()
        find_library(fftw3f_omp fftw3f_omp PATHS ${FFTW3F_LIBRARY_DIRS})
        if(fftw3f_omp)
            add_definitions(-DRT_FFTW3F_OMP)
            set(FFTW3F_LIBRARIES ${FFTW3F_LIBRARIES} ${fftw3f_omp})
        endif()
    endif()
endif()

# Find out whether we are building out of source:
get_filename_component(ABS_SOURCE_DIR "${PROJECT_SOURCE_DIR}" ABSOLUTE)
get_filename_component(ABS_BINARY_DIR "${CMAKE_BINARY_DIR}" ABSOLUTE)
set(OUT_OF_SOURCE_BUILD TRUE)
if(ABS_SOURCE_DIR STREQUAL ABS_BINARY_DIR)
    set(OUT_OF_SOURCE_BUILD FALSE)
    message(
        WARNING
            "You are performing an in-source build. This is discouraged. For an explanation and the advantages of out-of-source builds, please refer to http://www.cmake.org/Wiki/CMake_FAQ#What_is_an_.22out-of-source.22_build.3F"
    )
endif()

# Remove files which could require manual work. We will remove this after some
# time to have a clean build system without file modifications in the source
# tree again. (?)
set(OOSB_FILES
    "${PROJECT_SOURCE_DIR}/rtdata/rawtherapee.desktop"
    "${PROJECT_SOURCE_DIR}/rtgui/version.h"
    "${PROJECT_SOURCE_DIR}/rtgui/config.h"
    "${PROJECT_SOURCE_DIR}/AboutThisBuild.txt")
if(OUT_OF_SOURCE_BUILD)
    foreach(f ${OOSB_FILES})
        file(REMOVE "${f}")
    endforeach()
endif()

# Check for generated files in the source tree which should not be there when
# doing an out-of-source build. Without checking for this it might happen that
# old versions are used for the compilation:
if(OUT_OF_SOURCE_BUILD)
    foreach(f ${OOSB_FILES})
        if(EXISTS "${f}")
            message(
                SEND_ERROR
                    "Generated \"${f}\" found inside the source tree. Please remove it as it is a relic of the old build system and prevents valid compilation now."
            )
        endif()
    endforeach()
endif()

# Start generating AboutThisBuild.txt Set the platform bit-depth:
if(CMAKE_SIZEOF_VOID_P EQUAL 4)
    set(PROC_BIT_DEPTH 32 bits)
else()
    set(PROC_BIT_DEPTH 64 bits)
endif()

# Get compiler name and version.
get_filename_component(COMPILER_INFO ${CMAKE_C_COMPILER} NAME_WE)
set(COMPILER_INFO "${COMPILER_INFO} ${CMAKE_C_COMPILER_VERSION}")

# Get C++ and linker flags for rtengine (the GUI's C++ flags may have fewer
# flags):
set(CXX_FLAGS
    "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${UPPER_CMAKE_BUILD_TYPE}} ${RTENGINE_CXX_FLAGS}"
)
set(LFLAGS
    "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_${UPPER_CMAKE_BUILD_TYPE}}"
)

if(NOT APPLE)
    set(ABOUT_COMMAND_WITH_ARGS
        ${CMAKE_COMMAND} -DPROJECT_SOURCE_DIR:STRING=${PROJECT_SOURCE_DIR}
        -DCACHE_NAME_SUFFIX:STRING=${CACHE_NAME_SUFFIX}
        -DPROC_LABEL:STRING="${PROC_LABEL}"
        -DPROC_BIT_DEPTH:STRING="${PROC_BIT_DEPTH}"
        -DBUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
        -DGTKMM_VERSION:STRING=${GTKMM_VERSION}
        -DOPTION_OMP:STRING=${OPTION_OMP}
        -DWITH_MYFILE_MMAP:STRING=${WITH_MYFILE_MMAP}
        -DLENSFUN_VERSION:STRING=${LENSFUN_VERSION})
endif()

if(WIN32)
    list(
        APPEND
        ABOUT_COMMAND_WITH_ARGS
        -DSYSTEM:STRING=Windows
        -DCXX_FLAGS:STRING="${CXX_FLAGS}"
        -DLFLAGS:STRING="${LFLAGS}"
        -DCOMPILER_INFO:STRING="${COMPILER_INFO}"
        -DCMAKE_INSTALL_PREFIX:STRING="${CMAKE_INSTALL_PREFIX}"
        -DBIT_DEPTH:STRING="${CMAKE_SIZEOF_VOID_P}")
elseif(APPLE)
    set(
         ABOUT_COMMAND_WITH_ARGS
         cmake
             -DPROJECT_SOURCE_DIR:STRING=${PROJECT_SOURCE_DIR}
             -DCACHE_NAME_SUFFIX:STRING=${CACHE_NAME_SUFFIX}
             -P ${PROJECT_SOURCE_DIR}/UpdateInfo.cmake
             -DSYSTEM:STRING=Apple
             -DCXX_FLAGS:STRING=${CXX_FLAGS}
             -DLFLAGS:STRING=${LFLAGS}
             -DCOMPILER_INFO:STRING=${COMPILER_INFO})    
else()
    list(APPEND ABOUT_COMMAND_WITH_ARGS -DSYSTEM:STRING=Linux
         -DCXX_FLAGS:STRING=${CXX_FLAGS} -DLFLAGS:STRING=${LFLAGS}
         -DCOMPILER_INFO:STRING=${COMPILER_INFO})
endif()
if(NOT APPLE)
    list(APPEND
         ABOUT_COMMAND_WITH_ARGS
             -P "${PROJECT_SOURCE_DIR}/UpdateInfo.cmake")
endif()

add_custom_target(
    UpdateInfo ALL
    COMMAND ${ABOUT_COMMAND_WITH_ARGS}
    COMMENT "Creating AboutThisBuild.txt and other version-dependent files")

# End generating AboutThisBuild.txt

install(FILES AUTHORS.txt DESTINATION "${CREDITSDIR}")
install(FILES LICENSE.txt DESTINATION "${LICENCEDIR}")
install(FILES "${CMAKE_BINARY_DIR}/AboutThisBuild.txt"
        DESTINATION "${CREDITSDIR}")
install(
    FILES RELEASE_NOTES.txt
    DESTINATION "${CREDITSDIR}"
    OPTIONAL)

# The standard location for man pages in Linux is /usr/share/man Use "manpath"
# to see the search paths for man pages on your system.
if(BUILD_BUNDLE)
    if(APPLE)
        install(FILES "${PROJECT_SOURCE_DIR}/doc/manpage/rawtherapee.1"
                DESTINATION "${DATADIR}/man/man1")
    else()
        install(FILES "${PROJECT_SOURCE_DIR}/doc/manpage/rawtherapee.1"
                DESTINATION "${DATADIR}/share/man/man1")
    endif()
else()
    install(FILES "${PROJECT_SOURCE_DIR}/doc/manpage/rawtherapee.1"
            DESTINATION "${CMAKE_INSTALL_PREFIX}/share/man/man1")
endif()

if(WIN32)
    install(DIRECTORY "licenses" DESTINATION "${LICENCEDIR}")
endif()

if(UNIX)
    install(FILES com.rawtherapee.RawTherapee.appdata.xml
            DESTINATION "${APPDATADIR}")
endif()

# check whether the used version of lensfun has lfDatabase::LoadDirectory
include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_INCLUDES ${LENSFUN_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES)
foreach(l ${LENSFUN_LIBRARIES})
    if(LENSFUN_LIBRARY_DIRS)
        # the NO_DEFAULT_PATH is to make sure we find the lensfun version we
        # want, and not the system's one (e.g. if we have a custom version
        # installed in a non-standard location)
        find_library(_l ${l} PATHS ${LENSFUN_LIBRARY_DIRS} NO_DEFAULT_PATH)
    else()
        # LENSFUN_LIBRARY_DIRS can be empty if lensfun is installed in the
        # default path. In this case, adding NO_DEFAULT_PATH would make
        # find_library fail...
        find_library(_l ${l})
    endif()
    set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${_l})
endforeach()
check_cxx_source_compiles(
    "#include <lensfun.h>
int main()
{
    lfDatabase *db = 0;
    bool b = db->LoadDirectory(0);
    return 0;
}"
    LENSFUN_HAS_LOAD_DIRECTORY)

set(TCMALLOC_LIB_DIR)
if(ENABLE_TCMALLOC)
    if(TCMALLOC_LIB_DIR)
        find_library(TCMALLOC_LIBRARIES tcmalloc PATHS ${TCMALLOC_LIB_DIR}
                                                       NO_DEFAULT_PATH)
    else()
        find_library(TCMALLOC_LIBRARIES tcmalloc)
    endif()
    if(TCMALLOC_LIBRARIES)
        message(STATUS "using tcmalloc library in ${TCMALLOC_LIBRARIES}")
    else()
        set(TCMALLOC_LIBRARIES
            ""
            CACHE INTERNAL "" FORCE)
        message(STATUS "tcmalloc not found")
    endif()
else()
    set(TCMALLOC_LIBRARIES
        ""
        CACHE INTERNAL "" FORCE)
endif()

add_subdirectory(rtexif)
add_subdirectory(rtengine)
add_subdirectory(rtgui)
add_subdirectory(rtdata)
