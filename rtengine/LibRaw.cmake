# LibRaw has its own configuration script and uses make to build. Here we add
# the LibRaw configuration commands so they run during the CMake configuration.
# Also, add a target which always runs make.

set(LIBRAW_DIR "${CMAKE_CURRENT_BINARY_DIR}/libraw")
set(LIBRAW_LIB_DIR "${LIBRAW_DIR}/lib")
set(LIBRAW_PHANTOM_FILE "${LIBRAW_LIB_DIR}/phantom_file")
if(DEFINED ENV{SHELL})
    set(SHELL "$ENV{SHELL}")
else()
    set(SHELL "sh")
endif()

add_custom_target(
    LibRaw ALL
    DEPENDS ${LIBRAW_PHANTOM_FILE} # Ensures target always executes.
)

# Configuration flags.
set(CONFIGURE_FLAGS "--disable-examples")
set(LIBRAW_CXX_FLAGS "${CXX_FLAGS} -std=gnu++11 -Wno-error=unknown-pragmas")
# Let the configure script handle OpenMP flags.
string(REPLACE "${OpenMP_CXX_FLAGS}" "" LIBRAW_CXX_FLAGS "${LIBRAW_CXX_FLAGS}")
if(OPTION_OMP)
    set(CONFIGURE_FLAGS "${CONFIGURE_FLAGS} --enable-openmp")
else()
    set(CONFIGURE_FLAGS "${CONFIGURE_FLAGS} --disable-openmp")
endif()
set(CONFIGURE_FLAGS "${CONFIGURE_FLAGS} CC=\"${CMAKE_C_COMPILER}\"")
set(CONFIGURE_FLAGS "${CONFIGURE_FLAGS} CXX=\"${CMAKE_CXX_COMPILER}\"")
set(CONFIGURE_FLAGS "${CONFIGURE_FLAGS} CXXFLAGS=\"${LIBRAW_CXX_FLAGS}\"")

# Configuration commands.
message(STATUS "Configuring LibRaw")
execute_process(
    COMMAND cp -p -R "${CMAKE_CURRENT_SOURCE_DIR}/libraw" .
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    RESULT_VARIABLE PROCESS_RESULT
    COMMAND_ECHO STDOUT
)
if(PROCESS_RESULT AND NOT PROCESS_RESULT EQUAL 0)
    message(FATAL_ERROR "Could not copy LibRaw files into build directory")
endif()
execute_process(
    COMMAND "${SHELL}" -l -c "autoreconf -v --install"
    WORKING_DIRECTORY "${LIBRAW_DIR}"
    RESULT_VARIABLE PROCESS_RESULT
    COMMAND_ECHO STDOUT
)
if(PROCESS_RESULT AND NOT PROCESS_RESULT EQUAL 0)
    message(FATAL_ERROR "Could not generate LibRaw configuration script")
endif()
execute_process(
    COMMAND "${SHELL}" -l -c "./configure ${CONFIGURE_FLAGS}"
    WORKING_DIRECTORY "${LIBRAW_DIR}"
    RESULT_VARIABLE PROCESS_RESULT
    COMMAND_ECHO STDOUT
)
if(PROCESS_RESULT AND NOT PROCESS_RESULT EQUAL 0)
    execute_process(
        COMMAND cat config.log
        WORKING_DIRECTORY "${LIBRAW_DIR}"
        COMMAND_ECHO STDOUT
    )
    message(FATAL_ERROR "LibRaw configure failed")
endif()

# Build flags.
set(LIBRAW_MAKE_FLAGS "")
if(CMAKE_GENERATOR MATCHES ".*Makefiles.*")
    set(LIBRAW_MAKE_COMMAND "$(MAKE)")
else()
    # If not using Makefiles, set number of jobs equal to logical processors
    # count. Not necessary for make because of the jobserver.
    execute_process(
        COMMAND "${SHELL}" -l -c "nproc"
        OUTPUT_VARIABLE LOGICAL_PROCESSORS
        RESULT_VARIABLE PROCESS_RESULT
        ERROR_QUIET
    )
    if(PROCESS_RESULT AND NOT PROCESS_RESULT EQUAL 0)
        execute_process(
            COMMAND "${SHELL}" -l -c "sysctl -n hw.ncpu"
            OUTPUT_VARIABLE LOGICAL_PROCESSORS
            RESULT_VARIABLE PROCESS_RESULT
            ERROR_QUIET
        )
    endif()
    if(PROCESS_RESULT AND NOT PROCESS_RESULT EQUAL 0)
        execute_process(
            COMMAND "${SHELL}" -l -c "getconf _NPROCESSORS_ONLN"
            OUTPUT_VARIABLE LOGICAL_PROCESSORS
            RESULT_VARIABLE PROCESS_RESULT
            ERROR_QUIET
        )
    endif()
    if(PROCESS_RESULT AND NOT PROCESS_RESULT EQUAL 0)
        set(LOGICAL_PROCESSORS "1")
    endif()
    string(STRIP "${LOGICAL_PROCESSORS}" LOGICAL_PROCESSORS)
    set(LIBRAW_MAKE_FLAGS "${LIBRAW_MAKE_FLAGS} -j${LOGICAL_PROCESSORS}")

    set(LIBRAW_MAKE_COMMAND "make")
endif()

# Build commands.
add_custom_command(
    OUTPUT "${LIBRAW_PHANTOM_FILE}" "${LIBRAW_LIB_DIR}/.libs/libraw_r.a"
    COMMAND cp -p -R "${CMAKE_CURRENT_SOURCE_DIR}/libraw" ..
    COMMAND "${SHELL}" -l -c "${LIBRAW_MAKE_COMMAND} ${LIBRAW_MAKE_FLAGS}"
    COMMENT "Building LibRaw"
    WORKING_DIRECTORY libraw
    VERBATIM
)

# Add a `make clean-libraw` command because there's no good way to automatically
# clean the LibRaw build with `make`clean`.
add_custom_target(
    clean-libraw
    COMMAND make clean
    COMMAND rm -rf lib
    WORKING_DIRECTORY libraw
)

