IF (HALIDE_PATH)
    set(HALIDE_GENERATOR_TARGET "host")

    IF (OPTION_HALIDE_OPENCL)
        set(HALIDE_GENERATOR_TARGET ${HALIDE_GENERATOR_TARGET}-opencl)
    ENDIF (OPTION_HALIDE_OPENCL)
    link_directories("${HALIDE_PATH}/lib")
    include_directories("${HALIDE_PATH}/include")
    set(HALIDE_INCLUDE_DIRS "${HALIDE_PATH}/include")
    add_executable("igd_demosaic.generator" "${HALIDE_PATH}/tools/GenGen.cpp"
            "halide/generators/igd_demosaic_generator.cpp")
    target_compile_options("igd_demosaic.generator" PRIVATE "-std=c++11" "-fno-rtti")
    target_link_libraries("igd_demosaic.generator" "Halide" z)
ENDIF (HALIDE_PATH)

include_directories (${EXTRA_INCDIR} ${GTHREAD_INCLUDE_DIRS} ${GOBJECT_INCLUDE_DIRS} ${GLIB2_INCLUDE_DIRS}
    ${GLIBMM_INCLUDE_DIRS} ${IPTCDATA_INCLUDE_DIRS} ${LCMS_INCLUDE_DIRS} ${EXPAT_INCLUDE_DIRS} ${FFTW3F_INCLUDE_DIRS}
    ${GTKMM_INCLUDE_DIRS} ${GTK_INCLUDE_DIRS} ${HALIDE_INCLUDE_DIRS})
link_directories ("${PROJECT_SOURCE_DIR}/rtexif" ${EXTRA_LIBDIR} ${GTHREAD_LIBRARY_DIRS}
    ${GOBJECT_LIBRARY_DIRS} ${GLIB2_LIBRARY_DIRS} ${GLIBMM_LIBRARY_DIRS}
    ${IPTCDATA_LIBRARY_DIRS} ${LCMS_LIBRARY_DIRS} ${EXPAT_LIBRARY_DIRS} ${FFTW3F_LIBRARY_DIRS})

set (CAMCONSTSFILE "camconst.json")

set (RTENGINESOURCEFILES safegtk.cc colortemp.cc curves.cc flatcurves.cc diagonalcurves.cc dcraw.cc iccstore.cc color.cc
    dfmanager.cc ffmanager.cc rawimage.cc image8.cc image16.cc imagefloat.cc imagedata.cc imageio.cc improcfun.cc init.cc dcrop.cc
    loadinitial.cc procparams.cc rawimagesource.cc demosaic_algos.cc shmap.cc simpleprocess.cc refreshmap.cc
    fast_demo.cc amaze_demosaic_RT.cc CA_correct_RT.cc cfa_linedn_RT.cc green_equil_RT.cc hilite_recon.cc expo_before_b.cc
    stdimagesource.cc myfile.cc iccjpeg.cc hlmultipliers.cc improccoordinator.cc editbuffer.cc coord.cc
    processingjob.cc rtthumbnail.cc utils.cc labimage.cc slicer.cc cieimage.cc
    iplab2rgb.cc ipsharpen.cc iptransform.cc ipresize.cc ipvibrance.cc
    imagedimensions.cc jpeg_memsrc.cc jdatasrc.cc iimage.cc
    EdgePreservingDecomposition.cc cplx_wavelet_dec.cc FTblockDN.cc
    PF_correct_RT.cc previewimage.cc ipwavelet.cc
    dirpyr_equalizer.cc
    calc_distort.cc lcp.cc dcp.cc ipretinex.cc
    cJSON.c camconst.cc
    klt/convolve.cc klt/error.cc klt/klt.cc klt/klt_util.cc klt/pnmio.cc klt/pyramid.cc klt/selectGoodFeatures.cc
    klt/storeFeatures.cc klt/trackFeatures.cc klt/writeFeatures.cc
    clutstore.cc
    ciecam02.cc
    )

IF (HALIDE_PATH)
    add_definitions( -DHALIDE_ENABLED=1 )
    IF (OPTION_HALIDE_OPENCL)
        add_definitions( -DHALIDE_OPENCL=1 )
    ENDIF(OPTION_HALIDE_OPENCL)
    # Add the Halide code.
    set(RTENGINESOURCEFILES ${RTENGINESOURCEFILES} igd_halide_demosaic.cpp)
ENDIF (HALIDE_PATH)

include_directories (BEFORE "${CMAKE_CURRENT_BINARY_DIR}")

add_library (rtengine ${RTENGINESOURCEFILES})
add_dependencies (rtengine AboutFile)

IF (HALIDE_PATH)

    include("${CMAKE_CURRENT_SOURCE_DIR}/halide/HalideGenerator.cmake")

    halide_add_generator_dependency(TARGET rtengine
            GENERATOR_TARGET "igd_demosaic.generator"
            GENERATOR_NAME "igd_demosaic"
            GENERATED_FUNCTION "igd_demosaic"
            GENERATOR_ARGS "target=${HALIDE_GENERATOR_TARGET}")
    if (WIN32)
        # Windows options
    elseif(XCODE OR ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        # Mac OS X - Do nothing?
    else()
        # Unix
        set(HALIDE_LIBRARIES "dl pthread")
    endif()
ENDIF (HALIDE_PATH)

#It may be nice to store library version too
IF (BUILD_SHARED_LIBS)
    install (TARGETS rtengine DESTINATION ${LIBDIR})
ENDIF (BUILD_SHARED_LIBS)

set_target_properties (rtengine PROPERTIES COMPILE_FLAGS "${RTENGINE_CXX_FLAGS}")

target_link_libraries (rtengine rtexif ${EXTRA_LIB} ${GOBJECT_LIBRARIES} ${GTHREAD_LIBRARIES}
    ${GLIB2_LIBRARIES} ${GLIBMM_LIBRARIES} ${LCMS_LIBRARIES} ${EXPAT_LIBRARIES} ${FFTW3F_LIBRARIES} ${IPTCDATA_LIBRARIES}
    ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${TIFF_LIBRARIES} ${ZLIB_LIBRARIES} ${HALIDE_LIBRARIES})

install (FILES ${CAMCONSTSFILE} DESTINATION "${DATADIR}" PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
