add_library (rtexif rtexif.cc stdattribs.cc nikonattribs.cc canonattribs.cc
    pentaxattribs.cc fujiattribs.cc sonyminoltaattribs.cc olympusattribs.cc)

IF (WIN32)
    set_target_properties (rtexif PROPERTIES COMPILE_FLAGS " -ffast-math -fexpensive-optimizations")
ELSE (WIN32)
    set_target_properties (rtexif PROPERTIES COMPILE_FLAGS " -ffast-math -fexpensive-optimizations -fPIC")
ENDIF (WIN32)

IF (BUILD_SHARED_LIBS)
    INSTALL(TARGETS rtexif DESTINATION ${LIBDIR})
ENDIF (BUILD_SHARED_LIBS)
