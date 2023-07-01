add_library(rtexif STATIC
    canonattribs.cc
    fujiattribs.cc
    kodakattribs.cc
    nikonattribs.cc
    olympusattribs.cc
    panasonicattribs.cc
    pentaxattribs.cc
    rtexif.cc
    sonyminoltaattribs.cc
    stdattribs.cc
)

add_dependencies(rtexif UpdateInfo)

if(WIN32)
    include_directories(${EXTRA_INCDIR} ${GLIB2_INCLUDE_DIRS} ${GLIBMM_INCLUDE_DIRS} ${GTK_INCLUDE_DIRS} ${GTKMM_INCLUDE_DIRS} ${TIFF_INCLUDE_DIRS})
    link_directories(. "${PROJECT_SOURCE_DIR}/rtexif" ${EXTRA_LIBDIR} ${GLIB2_LIBRARY_DIRS} ${GLIBMM_LIBRARY_DIRS} ${GTK_LIBRARY_DIRS} ${GTKMM_LIBRARY_DIRS} ${LENSFUN_LIBRARY_DIRS} ${TIFF_LIBRARY_DIRS})
else()
    set_target_properties(rtexif PROPERTIES COMPILE_FLAGS " -fPIC")
    include_directories(${EXTRA_INCDIR} ${GLIB2_INCLUDE_DIRS} ${GLIBMM_INCLUDE_DIRS} ${GTK_INCLUDE_DIRS} ${GTKMM_INCLUDE_DIRS} ${TIFF_INCLUDE_DIRS})
    link_directories(${EXTRA_LIBDIR} ${GLIB2_LIBRARY_DIRS} ${GLIBMM_LIBRARY_DIRS} ${GTK_LIBRARY_DIRS} ${GTKMM_LIBRARY_DIRS} ${LENSFUN_LIBRARY_DIRS} ${TIFF_LIBRARY_DIRS})
endif()

include_directories(BEFORE "${CMAKE_CURRENT_BINARY_DIR}")

if(BUILD_SHARED_LIBS)
    install(TARGETS rtexif DESTINATION "${LIBDIR}")
endif()
