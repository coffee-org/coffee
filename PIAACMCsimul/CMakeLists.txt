# =====================================================================
# coffeePIAACMCsimul - PIAACMC Coronagraph Design & Simulation Suite
# =====================================================================

set(LIBNAME "coffeePIAACMCsimul")
set(SRCNAME "PIAACMCsimul")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(SOURCEFILES
    PIAACMCsimul.c
    init_piaacmcopticaldesign.c
    init_piaacmcopticalsystem.c
    init_piaacmcopticalsystem_elements.c
    PIAACMCsimul_conf_io.c
    PIAACMC_f_evalmask.c
    PIAACMCsimul_exec.c
    PIAACMCsimul_regularization.c
    PIAACMCsimul_linopt.c
    PIAACMCsimul_linopt_scan.c
    PIAACMCsimul_loadsavepiaacmcconf.c
    PIAACMCsimul_measure_transm_curve.c
    PIAACMCsimul_achromFPMsol_eval.c
    PIAACMCsimul_achromFPMsol_eval_zonezderivative.c
    exec_optimize_PIAA_shapes_fpmtransm.c
    PIAACMCsimul_CA2propCubeInt.c
    PIAACMCsimul_free.c
    PIAACMCsimul_computePSF.c
    PIAACMCsimul_computePSF_fast.c
    PIAACMCsimul_eval_poly_design.c
    PIAACMCsimul_eval_poly_sens.c
    PIAACMCsimul_run.c
    PIAACMCsimul_run_params.c
    exec_compute_image.c
    exec_computePSF_no_fpm.c
    FocalPlaneMask/exec_multizone_fpm_calib.c
    FocalPlaneMask/exec_multizone_fpm_calib_zones.c
    FocalPlaneMask/exec_optimize_fpm_zones.c
    FocalPlaneMask/exec_optimize_fpmtransmission.c
    FocalPlaneMask/FPM_process.c
    FocalPlaneMask/FPMresp_resample.c
    FocalPlaneMask/FPMresp_rmzones.c
    FocalPlaneMask/mkFocalPlaneMask.c
    FocalPlaneMask/mkFocalPlaneMask_cube.c
    FocalPlaneMask/mkFPM_zonemap.c
    FocalPlaneMask/rings2sectors.c
    LyotStop/exec_optimize_lyot_stop_position.c
    LyotStop/exec_optimize_lyot_stops_shapes_positions.c
    LyotStop/geomProp.c
    LyotStop/mkLyotMask.c
    LyotStop/mkSimpleLyotStop.c
    LyotStop/optimizeLyotStop.c
    LyotStop/optimizeLyotStop_eval.c
    PIAAshape/exec_optimize_PIAA_shapes.c
    PIAAshape/init_geomPIAA_rad.c
    PIAAshape/init_geomPIAA_rad_march.c
    PIAAshape/load2DRadialApodization.c
    PIAAshape/makePIAAshapes.c
    PIAAshape/mkPIAAMshapes_from_RadSag.c
    PIAAshape/PIAAshape_modes.c
    PIAAshape/setup_PIAA_shapes.c
)

set(INCLUDEFILES
    PIAACMCsimul.h
    exec_compute_image.h
    exec_computePSF_no_fpm.h
    exec_optimize_PIAA_shapes_fpmtransm.h
    init_piaacmcopticaldesign.h
    init_piaacmcopticalsystem.h
    PIAACMC_f_evalmask.h
    PIAACMCsimul_achromFPMsol_eval.h
    PIAACMCsimul_achromFPMsol_eval_zonezderivative.h
    PIAACMCsimul_CA2propCubeInt.h
    PIAACMCsimul_computePSF.h
    PIAACMCsimul_eval_poly_design.h
    PIAACMCsimul_exec.h
    PIAACMCsimul_loadsavepiaacmcconf.h
    PIAACMCsimul_measure_transm_curve.h
    PIAACMCsimul_run.h
    FocalPlaneMask/exec_multizone_fpm_calib.h
    FocalPlaneMask/exec_optimize_fpm_zones.h
    FocalPlaneMask/exec_optimize_fpmtransmission.h
    FocalPlaneMask/FPM_process.h
    FocalPlaneMask/FPMresp_resample.h
    FocalPlaneMask/FPMresp_rmzones.h
    FocalPlaneMask/mkFocalPlaneMask.h
    FocalPlaneMask/mkFPM_zonemap.h
    FocalPlaneMask/rings2sectors.h
    LyotStop/exec_optimize_lyot_stop_position.h
    LyotStop/exec_optimize_lyot_stops_shapes_positions.h
    LyotStop/geomProp.h
    LyotStop/mkLyotMask.h
    LyotStop/mkSimpleLyotStop.h
    LyotStop/optimizeLyotStop.h
    PIAAshape/exec_optimize_PIAA_shapes.h
    PIAAshape/init_geomPIAA_rad.h
    PIAAshape/load2DRadialApodization.h
    PIAAshape/makePIAAshapes.h
    PIAAshape/mkPIAAMshapes_from_RadSag.h
)

set(SCRIPTS
    scripts/coffee-optloop
    scripts/coffee-runclean
    scripts/coffee-runPIAACMC
    scripts/coffee-setconf
    scripts/coffee-syncscripts
    scripts/coffee-run
    scripts/coffee-runopt
    scripts/coffee-runPIAACMCdesign
    scripts/coffee-sim
    scripts/coffee-waitforfile1
    scripts/waitforfile
)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/scripts/coffee-piaacmc")
    list(APPEND SCRIPTS scripts/coffee-piaacmc)
endif()

project(lib_${LIBNAME}_project C)

find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED fftw3)
pkg_check_modules(FFTWF REQUIRED fftw3f)

if(COMMAND milk_add_runnable_module)
    milk_add_runnable_module(${LIBNAME}
        SOURCES ${SOURCEFILES}
        LINK_LIBS
            milkCOREMODiofits
            milkinfo
            milklinoptimtools
            milkimagebasic
            milkkdtree
            milkimagefilter
            milkimagegen
            milkfft
            milkstatistic
            milkZernikePolyn
            coffeeOptSystProp
            coffeecoronagraphs
            milkWFpropagate
            milkOpticsMaterials
            OpenMP::OpenMP_C
            GSL::gsl
            GSL::gslcblas
            ${FFTW_LIBRARIES}
            ${FFTWF_LIBRARIES}
    )
    milk_runnable_module_usage(${LIBNAME}
        INCLUDES
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/FocalPlaneMask>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/LyotStop>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/PIAAshape>
            $<INSTALL_INTERFACE:include>
            $<INSTALL_INTERFACE:include/${SRCNAME}>
    )
    milk_add_fpsexec_stubs(${LIBNAME})
else()
    include("${CMAKE_CURRENT_SOURCE_DIR}/../cmake/coffee_standalone.cmake")
    add_library(${LIBNAME} SHARED ${SOURCEFILES})
    coffee_configure_standalone_target(${LIBNAME}
        EXTRA_INCLUDES
            ${CMAKE_CURRENT_SOURCE_DIR}/FocalPlaneMask
            ${CMAKE_CURRENT_SOURCE_DIR}/LyotStop
            ${CMAKE_CURRENT_SOURCE_DIR}/PIAAshape
        MILK_LIBS
            milkinfo
            milklinoptimtools
            milkimagebasic
            milkkdtree
            milkimagefilter
            milkimagegen
            milkfft
            milkstatistic
            milkZernikePolyn
            milkWFpropagate
            milkOpticsMaterials
        LOCAL_TARGETS
            coffeeOptSystProp
            coffeecoronagraphs
    )
    install(TARGETS ${LIBNAME} DESTINATION lib)
endif()

install(FILES ${INCLUDEFILES} DESTINATION include/${SRCNAME})
install(PROGRAMS ${SCRIPTS} DESTINATION bin)
install(DIRECTORY examples/ DESTINATION share/coffee/PIAACMCsimul/examples)
