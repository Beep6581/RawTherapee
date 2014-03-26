add_library (rtexif rtexif.cc stdattribs.cc nikonattribs.cc canonattribs.cc
    pentaxattribs.cc fujiattribs.cc sonyminoltaattribs.cc olympusattribs.cc kodakattribs.cc)

IF (WIN32)
    set_target_properties (rtexif PROPERTIES COMPILE_FLAGS " -ffast-math -fexpensive-optimizations")
    include_directories (${EXTRA_INCDIR} ${GLIB2_INCLUDE_DIRS} ${GLIBMM_INCLUDE_DIRS}
        ${GTK_INCLUDE_DIRS} ${GTKMM_INCLUDE_DIRS})
    link_directories (. "${PROJECT_SOURCE_DIR}/rtexif" ${EXTRA_LIBDIR} ${GLIB2_LIBRARY_DIRS} ${GLIBMM_LIBRARY_DIRS}
        ${GTK_LIBRARY_DIRS} ${GTKMM_LIBRARY_DIRS})
    #set_target_properties (rth PROPERTIES LINK_FLAGS "-mwindows")
ELSE (WIN32)
    set_target_properties (rtexif PROPERTIES COMPILE_FLAGS " -ffast-math -fexpensive-optimizations -fPIC")
    include_directories (${EXTRA_INCDIR} ${GLIB2_INCLUDE_DIRS}
        ${GLIBMM_INCLUDE_DIRS} ${GTK_INCLUDE_DIRS} ${GTKMM_INCLUDE_DIRS})
    link_directories (${EXTRA_LIBDIR} ${GLIB2_LIBRARY_DIRS}
        ${GLIBMM_LIBRARY_DIRS} ${GTK_LIBRARY_DIRS} ${GTKMM_LIBRARY_DIRS})
ENDIF (WIN32)

include_directories (BEFORE "${CMAKE_CURRENT_BINARY_DIR}")

IF (BUILD_SHARED_LIBS)
    INSTALL(TARGETS rtexif DESTINATION ${LIBDIR})
ENDIF (BUILD_SHARED_LIBS)
