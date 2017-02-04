include_directories (${EXTRA_INCDIR} ${GTHREAD_INCLUDE_DIRS} ${GOBJECT_INCLUDE_DIRS} ${GLIB2_INCLUDE_DIRS}
    ${GLIBMM_INCLUDE_DIRS} ${IPTCDATA_INCLUDE_DIRS} ${LCMS_INCLUDE_DIRS} ${EXPAT_INCLUDE_DIRS} ${FFTW3F_INCLUDE_DIRS}
    ${GTKMM_INCLUDE_DIRS} ${GTK_INCLUDE_DIRS})
link_directories ("${PROJECT_SOURCE_DIR}/rtexif" ${EXTRA_LIBDIR} ${GTHREAD_LIBRARY_DIRS}
    ${GOBJECT_LIBRARY_DIRS} ${GLIB2_LIBRARY_DIRS} ${GLIBMM_LIBRARY_DIRS}
    ${IPTCDATA_LIBRARY_DIRS} ${LCMS_LIBRARY_DIRS} ${EXPAT_LIBRARY_DIRS} ${FFTW3F_LIBRARY_DIRS})
set (CAMCONSTSFILE "camconst.json")
set (RTENGINESOURCEFILES
	amaze_demosaic_RT.cc
	bilateral.cc
	CA_correct_RT.cc
	calc_distort.cc
	camconst.cc
	cfa_linedn_RT.cc
	ciecam02.cc
	cieimage.cc
	cJSON.c
	clutstore.cc
	color.cc
	colortemp.cc
	coord.cc
	cplx_wavelet_dec.cc
	curves.cc
	dcp.cc
	dcraw.cc
	dcrop.cc
	demosaic_algos.cc
	dfmanager.cc
	diagonalcurves.cc
	dirpyr_equalizer.cc
	EdgePreservingDecomposition.cc
	expo_before_b.cc
	fast_demo.cc
	ffmanager.cc
	flatcurves.cc
	FTblockDN.cc
	gauss.cc
	green_equil_RT.cc
	hilite_recon.cc
	iccjpeg.cc
	iccstore.cc
	iimage.cc
	image16.cc
	image8.cc
	imagedata.cc
	imagedimensions.cc
	imagefloat.cc
	imageio.cc
	improccoordinator.cc
	improcfun.cc
	init.cc
	iplab2rgb.cc
	ipresize.cc
	ipretinex.cc
	ipsharpen.cc
	iptransform.cc
	ipvibrance.cc
	ipwavelet.cc
	jdatasrc.cc
	jpeg_ijg/jpeg_memsrc.cc
	klt/convolve.cc
	klt/error.cc
	klt/klt.cc
	klt/klt_util.cc
	klt/pnmio.cc
	klt/pyramid.cc
	klt/selectGoodFeatures.cc
	klt/storeFeatures.cc
	klt/trackFeatures.cc
	klt/writeFeatures.cc
	labimage.cc
	lcp.cc
	loadinitial.cc
	myfile.cc
	PF_correct_RT.cc
	pipettebuffer.cc
	previewimage.cc
	processingjob.cc
	procparams.cc
	rawimage.cc
	rawimagesource.cc
	refreshmap.cc
	rtthumbnail.cc
	shmap.cc
	simpleprocess.cc
	slicer.cc
	stdimagesource.cc
	utils.cc
)
include_directories (BEFORE "${CMAKE_CURRENT_BINARY_DIR}")
add_library (rtengine ${RTENGINESOURCEFILES})
add_dependencies (rtengine UpdateInfo)
#It may be nice to store library version too
IF (BUILD_SHARED_LIBS)
    install (TARGETS rtengine DESTINATION ${LIBDIR})
ENDIF (BUILD_SHARED_LIBS)
set_target_properties (rtengine PROPERTIES COMPILE_FLAGS "${RTENGINE_CXX_FLAGS}")
target_link_libraries (rtengine rtexif ${EXTRA_LIB} ${GOBJECT_LIBRARIES} ${GTHREAD_LIBRARIES}
    ${GLIB2_LIBRARIES} ${GLIBMM_LIBRARIES} ${LCMS_LIBRARIES} ${EXPAT_LIBRARIES} ${FFTW3F_LIBRARIES} ${IPTCDATA_LIBRARIES}
    ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${TIFF_LIBRARIES} ${ZLIB_LIBRARIES})
install (FILES ${CAMCONSTSFILE} DESTINATION "${DATADIR}" PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)