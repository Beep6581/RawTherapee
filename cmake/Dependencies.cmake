include(FetchContent)

# Check for libraries:
find_package(PkgConfig)

# Wrap in function so updates to variables like CMAKE_CXX_FLAGS don't
# propagate out to other targets
function(rt_setup_dependencies)

    # fmt::fmt
    FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt
        GIT_TAG 11.1.4
        GIT_SHALLOW ON
    )

    # Add all FetchContent-declared libraries here.
    # Don't use FetchContent_Declare after this.
    FetchContent_MakeAvailable(
        fmt
    )

    # Gtk version shall be greater than 3.24.3 for fixed Hi-DPI support
    pkg_check_modules(GTK REQUIRED IMPORTED_TARGET gtk+-3.0>=3.24.3)
    pkg_check_modules(GTKMM REQUIRED IMPORTED_TARGET gtkmm-3.0>=3.24)

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

    pkg_check_modules(GLIB2 REQUIRED IMPORTED_TARGET glib-2.0>=2.48)
    pkg_check_modules(GLIBMM REQUIRED IMPORTED_TARGET glibmm-2.4>=2.48)
    pkg_check_modules(CAIROMM REQUIRED IMPORTED_TARGET cairomm-1.0)
    pkg_check_modules(GIO REQUIRED IMPORTED_TARGET gio-2.0>=2.48)
    pkg_check_modules(GIOMM REQUIRED IMPORTED_TARGET giomm-2.4>=2.48)
    pkg_check_modules(GTHREAD REQUIRED IMPORTED_TARGET gthread-2.0>=2.48)
    pkg_check_modules(GOBJECT REQUIRED IMPORTED_TARGET gobject-2.0>=2.48)
    pkg_check_modules(SIGC REQUIRED IMPORTED_TARGET sigc++-2.0>=2.3.1)

    pkg_check_modules(LENSFUN REQUIRED IMPORTED_TARGET lensfun>=0.2)
    pkg_check_modules(RSVG REQUIRED IMPORTED_TARGET librsvg-2.0>=2.52)

    pkg_check_modules(EXPAT REQUIRED IMPORTED_TARGET expat>=2.1)
    pkg_check_modules(FFTW3F REQUIRED IMPORTED_TARGET fftw3f)
    pkg_check_modules(IPTCDATA REQUIRED IMPORTED_TARGET libiptcdata)
    pkg_check_modules(LCMS REQUIRED IMPORTED_TARGET lcms2>=2.6)

    # By default little-cms2 uses 'register' keyword which is deprecated and
    # removed in c++17. This definition forces the library to not use it.
    target_compile_definitions(PkgConfig::LCMS INTERFACE CMS_NO_REGISTER_KEYWORD)

    find_package(ATOMIC)
    find_package(TIFF 4.0.4 REQUIRED)
    find_package(JPEG REQUIRED)
    find_package(PNG REQUIRED)
    find_package(ZLIB REQUIRED)

endfunction()
