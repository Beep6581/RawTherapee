#[=======================================================================[.rst:
FindTCMALLOC
------------

Find the tcmalloc library.

IMPORTED Targets
^^^^^^^^^^^^^^^^

This module provides the following :ref:`Imported Targets`:

``TCMALLOC::TCMALLOC``

  Target encapsulating the tcmalloc library usage requirements. It is available
  only when tcmalloc is found.

Result Variables
^^^^^^^^^^^^^^^^

This module sets the following variables:

``TCMALLOC_FOUND``
  True if tcmalloc is found.

``TCMALLOC_LIBRARIES``
  List of libraries needed to link to tcmalloc.

Hints
^^^^^

This module accepts the following variables:

``TCMALLOC_LIB_DIR``
  A user may set this variable to a tcmalloc installation root to help locate
  tcmalloc in custom installation paths.

``TCMALLOC_USE_MINIMAL``
  Set this variable to ``ON`` before calling ``find_package(TCMALLOC)`` to only
  link against tcmalloc_minimal.

Examples
^^^^^^^^

.. code-block:: cmake

  find_package(TCMALLOC)
  target_link_libraries(target PRIVATE TCMALLOC::TCMALLOC)
#]=======================================================================]

if(TCMALLOC_FIND_COMPONENTS AND NOT TCMALLOC_FIND_QUIETLY)
    message(AUTHOR_WARNING
        "TCMALLOC does not provide any COMPONENTS.  Calling\n"
        "  find_package(TCMALLOC COMPONENTS ...)\n"
        "will always fail."
    )
endif()

set(TCMALLOC_NAMES tcmalloc_minimal)
if(NOT TCMALLOC_USE_MINIMAL)
    list(APPEND TCMALLOC_NAMES tcmalloc)
endif()

# Allow TCMALLOC_LIBRARY to be set manually as the location of the library
if(NOT TCMALLOC_LIBRARY)
    if(TCMALLOC_LIB_DIR)
        find_library(
            TCMALLOC_LIBRARY ${TCMALLOC_NAMES}
            PATHS ${TCMALLOC_LIB_DIR}
            NO_DEFAULT_PATH
        )
    else()
        # On Debian libtcmalloc-minimal4 has no .so symlink, so let's adjust for that
        set(_tcmalloc_OLD_CMAKE_FIND_LIBRARY_SUFFIXES "${CMAKE_FIND_LIBRARY_SUFFIXES}")
        list(APPEND CMAKE_FIND_LIBRARY_SUFFIXES ".so.4")

        find_library(TCMALLOC_LIBRARY NAMES ${TCMALLOC_NAMES})

        set(CMAKE_FIND_LIBRARY_SUFFIXES "${_tcmalloc_OLD_CMAKE_FIND_LIBRARY_SUFFIXES}")
        unset(_tcmalloc_OLD_CMAKE_FIND_LIBRARY_SUFFIXES)
    endif()
endif()

unset(TCMALLOC_NAMES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TCMALLOC REQUIRED_VARS TCMALLOC_LIBRARY)

if(TCMALLOC_FOUND)
    if(NOT TCMALLOC_LIBRARIES)
        set(TCMALLOC_LIBRARIES ${TCMALLOC_LIBRARY})
    endif()

    if(NOT TARGET TCMALLOC::TCMALLOC)
        add_library(TCMALLOC::TCMALLOC UNKNOWN IMPORTED)
        set_target_properties(
            TCMALLOC::TCMALLOC
            PROPERTIES
            IMPORTED_LOCATION "${TCMALLOC_LIBRARY}"
        )
    endif()
endif()
