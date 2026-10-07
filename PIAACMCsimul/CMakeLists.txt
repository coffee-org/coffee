# library name
set(LIBNAME "coffeePIAACMCsimul")

set(SRCNAME "PIAACMCsimul")

message("")
message(" SRCNAME = ${SRCNAME} -> LIBNAME = ${LIBNAME}")


# list source files (.c) other than modulename.c
set(SOURCEFILES
	${SRCNAME}.c
    init_piaacmcopticaldesign.c
    init_piaacmcopticalsystem.c
    PIAACMC_f_evalmask.c
    PIAACMCsimul_exec.c
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
    PIAACMCsimul_run.c
    exec_compute_image.c
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
    init_piaacmcopticaldesign.h
    init_piaacmcopticalsystem.h
    PIAACMC_f_evalmask.h
    PIAACMCsimul_achromFPMsol_eval.h
    PIAACMCsimul_achromFPMsol_eval_zonezderivative.c
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
)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

project(lib_${LIBNAME}_project C)

find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED fftw3)
pkg_check_modules(FFTWF REQUIRED fftw3f)

add_library(${LIBNAME} SHARED ${SOURCEFILES})

set_target_properties(${LIBNAME} PROPERTIES
    LINKER_LANGUAGE C
    INTERPROCEDURAL_OPTIMIZATION OFF
)

if(TARGET CLIcore)
    # In-tree build under milk framework
    target_include_directories(${LIBNAME} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/FocalPlaneMask
        ${CMAKE_CURRENT_SOURCE_DIR}/LyotStop
        ${CMAKE_CURRENT_SOURCE_DIR}/PIAAshape
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/coremods
        ${CMAKE_SOURCE_DIR}/plugins/milk-extra-src
        $<TARGET_PROPERTY:CLIcore,INTERFACE_INCLUDE_DIRECTORIES>
    )

    set(OPT_WF_LIBS milkWFpropagate milkOpticsMaterials)

    target_link_libraries(${LIBNAME} PUBLIC
        CLIcore
        ImageStreamIO
        milkCOREMODmemory
        milkCOREMODarith
        milkCOREMODiofits
        milkCOREMODtools
        milkinfo
        milklinoptimtools
        milkimagebasic
        milkimagefilter
        milkimagegen
        milkfft
        milkstatistic
        milkZernikePolyn
        coffeeOptSystProp
        coffeecoronagraphs
        ${OPT_WF_LIBS}
        OpenMP::OpenMP_C
        GSL::gsl
        GSL::gslcblas
        ${FFTW_LIBRARIES}
        ${FFTWF_LIBRARIES}
        m
    )
else()
    # Out-of-tree / standalone build
    set(MILK_SOURCE_CANDIDATES
        ${MILK_SOURCE_DIR}
        $ENV{MILK_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/../../milk-framework-dev
        $ENV{MILK_ROOT}
        ${CMAKE_CURRENT_SOURCE_DIR}/../../milk
        /usr/local/milk
    )
    set(MILK_SOURCE_FOUND "")
    foreach(cand ${MILK_SOURCE_CANDIDATES})
        if(cand AND EXISTS "${cand}/src/coremods")
            get_filename_component(MILK_SOURCE_FOUND "${cand}" ABSOLUTE)
            break()
        endif()
    endforeach()

    find_path(MILK_INCLUDE_DIR
        NAMES CommandLineInterface/CLIcore.h CLIcore.h
        HINTS
            ${MILK_ROOT}/include
            $ENV{MILK_ROOT}/include
            /usr/local/milk/include
            /usr/local/include
    )

    if(NOT DEFINED MILK_STANDALONE_INCLUDES)
        set(MILK_STANDALONE_INCLUDES
            ${CMAKE_CURRENT_SOURCE_DIR}/..
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${CMAKE_CURRENT_SOURCE_DIR}/FocalPlaneMask
            ${CMAKE_CURRENT_SOURCE_DIR}/LyotStop
            ${CMAKE_CURRENT_SOURCE_DIR}/PIAAshape
            /usr/local/milk/include
            /usr/local/include
        )

        if(MILK_SOURCE_FOUND)
            list(APPEND MILK_STANDALONE_INCLUDES
                ${MILK_SOURCE_FOUND}/src
                ${MILK_SOURCE_FOUND}/src/coremods
                ${MILK_SOURCE_FOUND}/plugins/milk-extra-src
                ${CMAKE_CURRENT_SOURCE_DIR}/../../OpticsMaterials
                ${CMAKE_CURRENT_SOURCE_DIR}/../../WFpropagate
                ${MILK_SOURCE_FOUND}/plugins/OpticsMaterials
                ${MILK_SOURCE_FOUND}/plugins/WFpropagate
                ${MILK_SOURCE_FOUND}/plugins/milk-extra-src/OpticsMaterials
                ${MILK_SOURCE_FOUND}/plugins/milk-extra-src/WFpropagate
                ${MILK_SOURCE_FOUND}/src/cli
                ${MILK_SOURCE_FOUND}/src/cli/CLIcore
                ${MILK_SOURCE_FOUND}/src/cli/libmilkscript
                ${MILK_SOURCE_FOUND}/src/engine
                ${MILK_SOURCE_FOUND}/src/engine/libfps
                ${MILK_SOURCE_FOUND}/src/engine/libmilkcommon
                ${MILK_SOURCE_FOUND}/src/engine/libprocessinfo
                ${MILK_SOURCE_FOUND}/src/engine/libmilkdata
                ${MILK_SOURCE_FOUND}/src/engine/ImageStreamIO
            )
            if(EXISTS "${MILK_SOURCE_FOUND}/_build/config")
                list(APPEND MILK_STANDALONE_INCLUDES "${MILK_SOURCE_FOUND}/_build/config")
            endif()
        endif()

        if(MILK_INCLUDE_DIR)
            list(APPEND MILK_STANDALONE_INCLUDES
                ${MILK_INCLUDE_DIR}
                ${MILK_INCLUDE_DIR}/CommandLineInterface
                ${MILK_INCLUDE_DIR}/libmilkcommon
                ${MILK_INCLUDE_DIR}/libfps
                ${MILK_INCLUDE_DIR}/libmilkdata
                ${MILK_INCLUDE_DIR}/ImageStreamIO
            )
        endif()
    endif()

    target_include_directories(${LIBNAME} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/FocalPlaneMask
        ${CMAKE_CURRENT_SOURCE_DIR}/LyotStop
        ${CMAKE_CURRENT_SOURCE_DIR}/PIAAshape
        ${MILK_STANDALONE_INCLUDES}
    )

    set(MILK_LIB_NAMES
        CLIcore
        ImageStreamIO
        milkCOREMODmemory
        milkCOREMODarith
        milkCOREMODiofits
        milkCOREMODtools
        milkinfo
        milklinoptimtools
        milkimagebasic
        milkimagefilter
        milkimagegen
        milkfft
        milkstatistic
        milkZernikePolyn
        milkWFpropagate
        milkOpticsMaterials
    )

    set(MILK_SEARCH_HINTS
        ${CMAKE_CURRENT_SOURCE_DIR}/../../OpticsMaterials/_build
        ${CMAKE_CURRENT_SOURCE_DIR}/../../WFpropagate/_build
    )
    if(MILK_SOURCE_FOUND)
        list(APPEND MILK_SEARCH_HINTS
            ${MILK_SOURCE_FOUND}/_build/lib
            ${MILK_SOURCE_FOUND}/_build/plugins/OpticsMaterials
            ${MILK_SOURCE_FOUND}/_build/plugins/WFpropagate
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/OpticsMaterials
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/WFpropagate
            ${MILK_SOURCE_FOUND}/_build/src/cli/CLIcore
            ${MILK_SOURCE_FOUND}/_build/src/engine/ImageStreamIO
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/fft
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/image_gen
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/image_basic
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/image_filter
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/linopt_imtools
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/info
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/statistic
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/ZernikePolyn
        )
    endif()
    list(APPEND MILK_SEARCH_HINTS
        ${MILK_ROOT}/lib
        $ENV{MILK_ROOT}/lib
        /usr/local/milk/lib
        /usr/local/lib
    )

    set(MILK_LIBRARIES "")
    foreach(lib ${MILK_LIB_NAMES})
        find_library(MILK_${lib}_LIB NAMES ${lib} HINTS ${MILK_SEARCH_HINTS})
        if(MILK_${lib}_LIB)
            list(APPEND MILK_LIBRARIES ${MILK_${lib}_LIB})
        else()
            list(APPEND MILK_LIBRARIES ${lib})
        endif()
    endforeach()

    target_link_libraries(${LIBNAME} PUBLIC
        ${MILK_LIBRARIES}
        coffeeOptSystProp
        coffeecoronagraphs
        OpenMP::OpenMP_C
        GSL::gsl
        GSL::gslcblas
        ${FFTW_LIBRARIES}
        ${FFTWF_LIBRARIES}
        m
    )
endif()

install(TARGETS ${LIBNAME} EXPORT milkTargets DESTINATION lib)
install(FILES ${INCLUDEFILES} DESTINATION include/${SRCNAME})
install(PROGRAMS ${SCRIPTS} DESTINATION bin)


