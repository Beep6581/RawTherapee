# - Try to find KLT
# Once done this will define
#
#  KLT_FOUND - system has KLT
#  KLT_INCLUDE_DIRS - the KLT include directory
#  KLT_LIBRARIES - Link these to use KLT
#  KLT_DEFINITIONS - Compiler switches required for using KLT
#
#  Copyright (c) 2009 Andreas Schneider <mail@cynapses.org>
#  updated for KLT by Dan Horák
#
#  Redistribution and use is allowed according to the terms of the New
#  BSD license.
#  For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#


if (KLT_LIBRARIES AND KLT_INCLUDE_DIRS)
  # in cache already
  set(KLT_FOUND TRUE)
else (KLT_LIBRARIES AND KLT_INCLUDE_DIRS)

  find_path(KLT_INCLUDE_DIR
    NAMES
        klt.h
    PATHS
        /usr/include
        /usr/local/include
        /opt/local/include
        /sw/include
        ${CMAKE_INSTALL_PREFIX}/include
    PATH_SUFFIXES
        klt
  )
  mark_as_advanced(KLT_INCLUDE_DIR)

  find_library(KLT_LIBRARY
    NAMES
        klt
    PATHS
        /usr/lib64
        /usr/lib
        /usr/local/lib
        /opt/local/lib
        /sw/lib
        ${CMAKE_INSTALL_PREFIX}/lib
  )
  mark_as_advanced(KLT_LIBRARY)

  set(KLT_INCLUDE_DIRS
    ${KLT_INCLUDE_DIR}
  )

  set(KLT_LIBRARIES
    ${KLT_LIBRARY}
  )

  if (KLT_INCLUDE_DIRS AND KLT_LIBRARIES)
     set(KLT_FOUND TRUE)
  endif (KLT_INCLUDE_DIRS AND KLT_LIBRARIES)

  if (KLT_FOUND)
    if (NOT KLT_FIND_QUIETLY)
      message(STATUS "Found KLT: ${KLT_LIBRARIES}")
    endif (NOT KLT_FIND_QUIETLY)
  else (KLT_FOUND)
      message(STATUS "KLT not found.")
  endif (KLT_FOUND)

  mark_as_advanced(KLT_INCLUDE_DIRS KLT_LIBRARIES)

endif (KLT_LIBRARIES AND KLT_INCLUDE_DIRS)

