# Use the 'Debug'          build type to have a non optimized,        with    debugging information, with a console executable
# Use the 'Release'        build type to have an    optimized,        without debugging information, console free executable
# Use the 'RelWithDebInfo' build type to have an    optimized,        without debugging information, with a console executable
# Use the 'MinSizeRel'     build type to have the smallest possible,  without debugging information, console free executable
#set(CMAKE_BUILD_TYPE Release CACHE STRING "Between: None Debug Release RelWithDebInfo MinSizeRel.")

set(CMAKE_INSTALL_PREFIX ./${CMAKE_BUILD_TYPE} CACHE PATH "Libraries installation path")
set(DATADIR . CACHE PATH "Data installation path")
set(BINDIR . CACHE PATH "Binaries installation path")
set(LIBDIR . CACHE PATH "Libraries installation path")
set(DOCDIR ./doc CACHE PATH "Documentation installation path")
set(CREDITSDIR . CACHE PATH "Credit file installation path")
set(LICENCEDIR . CACHE PATH "Licence file installation path")

set(BUILD_SHARED OFF CACHE BOOL "Should RT generate shared libraries")
set(OPTION_OMP ON CACHE BOOL "Use OpenMP to speedup the preview and batch processing")
# set WITH_MYFILE_MMAP to OFF if you experience crash with thumbnail creation (it should be slower, but more reliable)
set(WITH_MYFILE_MMAP ON CACHE BOOL "Use the MMAP mechanism to speedup thumbnail creations")

set(CACHE_NAME_SUFFIX "" CACHE STRING "RawTherapee's cache folder suffix (leave empty to use the default suffix, i.e. latesttag)")

# This line will let you chose the target number, and the associated processor
set(PROC_TARGET_NUMBER 0 CACHE STRING "Target Processor")

# If you want to force the target processor name when PROC_TARGET_NUMBER = 0 or 2,
# uncomment the next line and replace labelWithoutQuotes by its value
#set(PROC_LABEL labelWithoutQuotes CACHE STRING "Target Processor label")

# Important: MinGW-w64 user may need to specify the -m32 or -m64 flag in CMAKE_CXX_FLAGS,
# CMAKE_C_FLAGS and CMAKE_EXE_LINKER_FLAGS to select between 32/64bit build
set(CMAKE_CXX_FLAGS "-mwin32 -m64 -mthreads -msse2" CACHE STRING "Compiler options for C++ source files")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g2" CACHE STRING "Compiler options for C++ source files and Debug target")
set(CMAKE_CXX_FLAGS_RELEASE "-Wno-aggressive-loop-optimizations -DNDEBUG -O3" CACHE STRING "Compiler options for C++ source files and Release target")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-DNDEBUG -Os" CACHE STRING "Compiler options for C++ source files and MinSizeRel target")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g" CACHE STRING "Compiler options for C++ source files and RelWithDebInfo target")

# Uncomment the next line and set the right value to override the default value (special compiling flags for RTEngine)
#set(RTENGINE_CXX_FLAGS "-funroll-loops" CACHE STRING "Special compilation flags for RTEngine")

set(CMAKE_C_FLAGS "-mwin32 -m64 -mthreads -msse2" CACHE STRING "Compiler options for C source files")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g2" CACHE STRING "Compiler options for C source files and Debug target")
set(CMAKE_C_FLAGS_RELEASE "-DNDEBUG -O2" CACHE STRING "Compiler options for C source files and Release target")
set(CMAKE_C_FLAGS_MINSIZEREL "-DNDEBUG -Os" CACHE STRING "Compiler options for C source files and MinSizeRel target")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g" CACHE STRING "Compiler options for C source files and RelWithDebInfo target")

set(CMAKE_EXE_LINKER_FLAGS "-m64 -mthreads -static-libgcc" CACHE STRING "Linker options")
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-O0" CACHE STRING "Linkage options for the Debug target")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-s -O3" CACHE STRING "Linkage options for the Release target")
set(CMAKE_EXE_LINKER_FLAGS_MINSIZEREL "-s -Os" CACHE STRING "Linkage options for the MinSizeRel target")
set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "-s -O2" CACHE STRING "Linkage options for the RelWithDebInfo target")
