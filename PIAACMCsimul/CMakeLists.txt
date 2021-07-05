# library name
set(LIBNAME "coffeePIAACMCsimul")

set(SRCNAME "PIAACMCsimul")

message("")
message(" SRCNAME = ${SRCNAME} -> LIBNAME = ${LIBNAME}")


# list source files (.c) other than modulename.c
set(SOURCEFILES
	${SRCNAME}.c
    PIAACMC_f_evalmask.c
    PIAACMCsimul_loadsavepiaacmcconf.c
    PIAACMCsimul_measure_transm_curve.c
    PIAACMCsimul_achromFPMsol_eval.c
    PIAACMCsimul_achromFPMsol_eval_zonezderivative.c
    PIAACMCsimul.c
    exec_optimize_PIAA_shapes_fpmtransm.c
    PIAACMCsimul_CA2propCubeInt.c
    PIAACMCsimul_free.c
    PIAACMCsimul_computePSF.c
    PIAACMCsimul_eval_poly_design.c
    PIAACMCsimul_init.c
    PIAACMCsimul_exec.c
    PIAACMCsimul_run.c
    exec_compute_image.c
    PIAACMCsimul_initpiaacmcconf.c
    exec_computePSF_no_fpm.c
    FocalPlaneMask/exec_multizone_fpm_calib.c
    FocalPlaneMask/exec_optimize_fpm_zones.c
    FocalPlaneMask/exec_optimize_fpmtransmission.c
    FocalPlaneMask/FPM_process.c
    FocalPlaneMask/FPMresp_resample.c
    FocalPlaneMask/FPMresp_rmzones.c
    FocalPlaneMask/mkFocalPlaneMask.c
    FocalPlaneMask/mkFPM_zonemap.c
    FocalPlaneMask/rings2sectors.c
    LyotStop/exec_optimize_lyot_stop_position.c
    LyotStop/exec_optimize_lyot_stops_shapes_positions.c
    LyotStop/geomProp.c
    LyotStop/mkLyotMask.c
    LyotStop/mkSimpleLyotStop.c
    LyotStop/optimizeLyotStop.c
    PIAAshape/exec_optimize_PIAA_shapes.c
    PIAAshape/init_geomPIAA_rad.c
    PIAAshape/load2DRadialApodization.c
    PIAAshape/makePIAAshapes.c
    PIAAshape/mkPIAAMshapes_from_RadSag.c
)

set(INCLUDEFILES
	${SRCNAME}.h
    exec_compute_image.h
    exec_computePSF_no_fpm.h
    exec_optimize_PIAA_shapes_fpmtransm.h
    PIAACMC_f_evalmask.h
    PIAACMCsimul_achromFPMsol_eval.h
    PIAACMCsimul_achromFPMsol_eval_zonezderivative.c
    PIAACMCsimul_CA2propCubeInt.h
    PIAACMCsimul_computePSF.h
    PIAACMCsimul_eval_poly_design.h
    PIAACMCsimul_exec.h
    PIAACMCsimul_init.h
    PIAACMCsimul_initpiaacmcconf.h
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
)

set(LINKLIBS
	CLIcore
    milkinfo
    milkOpticsMaterials
    milklinoptimtools
    milkimagebasic
    milkimagefilter
)

# DEFAULT SETTINGS
# Do not change unless needed
# =====================================================================

project(lib_${LIBNAME}_project)

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${PROJECT_SOURCE_DIR}/..")


# Library can be compiled from multiple source files
# Convention: the main souce file is named <libname>.c
#
add_library(${LIBNAME} SHARED ${SOURCEFILES})

target_include_directories(${LIBNAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(${LIBNAME} PRIVATE ${LINKLIBS})

install(TARGETS ${LIBNAME} DESTINATION lib)
install(FILES ${INCLUDEFILES} DESTINATION include/${SRCNAME})

install(PROGRAMS ${SCRIPTS} DESTINATION bin)
