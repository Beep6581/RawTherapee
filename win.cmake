# Use the 'Debug'          build type to have a non optimized,        with    debugging information, with a console executable.
# Use the 'Release'        build type to have an    optimized,        without debugging information, console-free executable.
# Use the 'RelWithDebInfo' build type to have an    optimized,        without debugging information, with a console executable.
# Use the 'MinSizeRel'     build type to have the smallest possible,  without debugging information, console-free executable.
#
#set(CMAKE_BUILD_TYPE Release CACHE STRING "Between: None Debug Release RelWithDebInfo MinSizeRel.")

set(CMAKE_INSTALL_PREFIX ./Builds/${CMAKE_BUILD_TYPE} CACHE PATH "Library installation path")
set(DATADIR . CACHE PATH "Data installation path")
set(BINDIR . CACHE PATH "Binaries installation path")
set(LIBDIR . CACHE PATH "Library installation path")
set(DOCDIR ./doc CACHE PATH "Documentation installation path")
set(CREDITSDIR . CACHE PATH "Credit file installation path")
set(LICENCEDIR . CACHE PATH "Licence file installation path")

set(BUILD_SHARED OFF CACHE BOOL "Should shared libraries be generated")
set(OPTION_OMP ON CACHE BOOL "Use OpenMP to speedup the preview and batch processing")
# set WITH_MYFILE_MMAP to OFF if you experience crashes caused by thumbnail creation (it will be slower but more reliable).
set(WITH_MYFILE_MMAP ON CACHE BOOL "Use the MMAP mechanism to speedup thumbnail creation")

set(CACHE_NAME_SUFFIX "" CACHE STRING "RawTherapee's cache folder suffix (leave empty to use default from UpdateIndo.cmake)")

# Choose the target processor's number from ProcessorTargets.cmake.
set(PROC_TARGET_NUMBER 0 CACHE STRING "Target processor")

# If you want to force the target processor name when PROC_TARGET_NUMBER = 0 or 2,
# uncomment the next line and replace labelWithoutQuotes with a string.
#set (PROC_LABEL labelWithoutQuotes CACHE STRING "Target processor label")

# Important: MinGW-w64 users may need to specify the -m32 or -m64 flag in CMAKE_CXX_FLAGS,
# CMAKE_C_FLAGS and CMAKE_EXE_LINKER_FLAGS to select between 32-/64-bit builds.
set(CMAKE_CXX_FLAGS "-mwin32 -mthreads" CACHE STRING "Compiler options for C++ source files")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g2" CACHE STRING "Compiler options for C++ source files and Debug target")
set(CMAKE_CXX_FLAGS_RELEASE "-mwindows -DNDEBUG -O2" CACHE STRING "Compiler options for C++ source files and Release target")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-mwindows -DNDEBUG -Os" CACHE STRING "Compiler options for C++ source files and MinSizeRel target")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g" CACHE STRING "Compiler options for C++ source files and RelWithDebInfo target")

set(CMAKE_C_FLAGS "-mwin32 -mthreads" CACHE STRING "Compiler options for C source files")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g2" CACHE STRING "Compiler options for C source files and Debug target")
set(CMAKE_C_FLAGS_RELEASE "-mwindows -DNDEBUG -O2" CACHE STRING "Compiler options for C source files and Release target")
set(CMAKE_C_FLAGS_MINSIZEREL "-mwindows -DNDEBUG -Os" CACHE STRING "Compiler options for C source files and MinSizeRel target")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g" CACHE STRING "Compiler options for C source files and RelWithDebInfo target")

set(CMAKE_EXE_LINKER_FLAGS "-mwin32 -mthreads -static-libgcc -Wl,--large-address-aware,--verbose" CACHE STRING "Linker options")
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-O0" CACHE STRING "Linkage options for the Debug target")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-mwindows -s -O2" CACHE STRING "Linkage options for the Release target")
set(CMAKE_EXE_LINKER_FLAGS_MINSIZEREL "-mwindows -s -Os" CACHE STRING "Linkage options for the MinSizeRel target")
set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "-s -O2" CACHE STRING "Linkage options for the RelWithDebInfo target")
