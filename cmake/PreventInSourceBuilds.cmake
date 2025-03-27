function(rt_force_out_of_source_builds)
  get_filename_component(SRC_DIR "${CMAKE_SOURCE_DIR}" REALPATH)
  get_filename_component(BIN_DIR "${CMAKE_BINARY_DIR}" REALPATH)

  if("${SRC_DIR}" STREQUAL "${BIN_DIR}")
    message("############################################################")
    message("Warning: In-source builds are disabled.")
    message("Please create and run cmake from a separate build directory.")
    message("############################################################")
    message(FATAL_ERROR "Quitting configuration")
  endif()
endfunction()

rt_force_out_of_source_builds()
