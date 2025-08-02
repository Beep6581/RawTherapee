function(rt_force_out_of_source_builds)
    get_filename_component(SRC_DIR "${CMAKE_SOURCE_DIR}" REALPATH)
    get_filename_component(BIN_DIR "${CMAKE_BINARY_DIR}" REALPATH)

    if("${SRC_DIR}" STREQUAL "${BIN_DIR}")
        message(WARNING "############################################################")
        message(WARNING "Warning: In-source builds are disabled.")
        message(WARNING "Please create and run cmake from a separate build directory.")
        message(WARNING "############################################################")
        message(FATAL_ERROR "Quitting configuration...")
    endif()
endfunction()

macro(rt_basic_build_settings)
    # Set a default build type if none was specified
    if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        message(STATUS "Setting build type to 'Release' as none was specified.")
        set(CMAKE_BUILD_TYPE
            Release
            CACHE STRING "Choose the type of build." FORCE
        )

        # Set the possible values of build type for cmake-gui, ccmake
        set_property(
            CACHE CMAKE_BUILD_TYPE
            PROPERTY STRINGS
            "Debug"
            "Release"
            "MinSizeRel"
            "RelWithDebInfo"
        )
    endif()

    if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
        set(IS_DEBUG ON)
    else()
        set(IS_DEBUG OFF)
    endif()

    # Generate compile_commands.json to make it easier to work with clang based tools
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endmacro()
