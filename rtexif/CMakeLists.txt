include_directories (.)
add_library (rtexif rtexif.cc stdattribs.cc nikonattribs.cc canonattribs.cc
    pentaxattribs.cc fujiattribs.cc sonyminoltaattribs.cc olympusattribs.cc)

IF (WIN32)
    set_target_properties (rtexif PROPERTIES COMPILE_FLAGS "-O3 -ffast-math -fexpensive-optimizations")
ELSE (WIN32)
    set_target_properties (rtexif PROPERTIES COMPILE_FLAGS "-O3 -ffast-math -fexpensive-optimizations -fPIC")
ENDIF (WIN32)
 