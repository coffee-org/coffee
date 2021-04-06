set(LIBNAME "coffeePIAACMCsimul")
set(SRCNAME "PIAACMCsimul")


set(SOURCEFILES
	${SRCNAME}.c
    PIAACMC_f_evalmask.c
    PIAACMCsimul_exec_multizone_fpm_calib.c
    PIAACMCsimul_loadsavepiaacmcconf.c
    PIAACMC_FPM_process.c
    PIAACMCsimul_exec_optimize_fpmtransmission.c
    PIAACMCsimul_makePIAAshapes.c
    PIAACMC_FPMresp_resample.c
    PIAACMCsimul_exec_optimize_fpm_zones.c
    PIAACMCsimul_measure_transm_curve.c
    PIAACMC_FPMresp_rmzones.c
    PIAACMCsimul_exec_optimize_lyot_stop_position.c
    PIAACMCsimul_mkFocalPlaneMask.c
    PIAACMCsimul_achromFPMsol_eval.c
    PIAACMCsimul_exec_optimize_lyot_stops_shapes_positions.c
    PIAACMCsimul_mkFPM_zonemap.c
    PIAACMCsimul_achromFPMsol_eval_zonezderivative.c
    PIAACMCsimul_exec_optimize_PIAA_shapes.c
    PIAACMCsimul_mkLyotMask.c
    PIAACMCsimul.c
    PIAACMCsimul_exec_optimize_PIAA_shapes_fpmtransm.c
    PIAACMCsimul_mkPIAAMshapes_from_RadSag.c
    PIAACMCsimul_CA2propCubeInt.c
    PIAACMCsimul_free.c
    PIAACMCsimul_mkSimpleLyotStop.c
    PIAACMCsimul_computePSF.c
    PIAACMCsimul_geomProp.c
    PIAACMCsimul_optimizeLyotStop.c
    PIAACMCsimul_eval_poly_design.c
    PIAACMCsimul_init.c
    PIAACMCsimul_rings2sectors.c
    PIAACMCsimul_exec.c
    PIAACMCsimul_init_geomPIAA_rad.c
    PIAACMCsimul_run.c
    PIAACMCsimul_exec_compute_image.c
    PIAACMCsimul_initpiaacmcconf.c
    PIAACMCsimul_exec_computePSF_no_fpm.c
    PIAACMCsimul_load2DRadialApodization.c
)

set(INCLUDEFILES
	${SRCNAME}.h
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
