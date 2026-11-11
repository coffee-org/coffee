# =====================================================================
# coffeecoronagraphs - Coronagraph Design & Propagation Algorithms
# =====================================================================

cmake_minimum_required(VERSION 3.16)

set(LIBNAME "coffeecoronagraphs")
set(SRCNAME "coronagraphs")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(SOURCEFILES
    coronagraphs.c
    coronagraphs_prolate.c
    coronagraphs_prolate_dft.c
    coronagraphs_prolate_update.c
    coronagraphs_apofit_eval.c
    coronagraphs_apofit.c
    coronagraphs_aplc_compile.c
    coronagraphs_piaa_init.c
    coronagraphs_piaa_remap.c
    coronagraphs_piaa_beam.c
    coronagraphs_masks_init.c
    coronagraphs_sim_phase.c
    coronagraphs_sim_bandlim.c
    coronagraphs_sim_vortex_rrpm.c
    coronagraphs_sim_apodized.c
    coronagraphs_sim_piaac.c
    coronagraphs_sim_aplc.c
    coronagraphs_sim_aplc_prop.c
    coronagraphs_eval.c
    coronagraphs_transm.c
    coronagraphs_limitcoeff.c
    coronagraphs_scan_compute.c
    coronagraphs_scan.c
    coronagraphs_cli.c
)

set(INCLUDEFILES
    coronagraphs.h
)

set(SCRIPTS
    scripts/coffee
)

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
            milkimagegen
            milkfft
            milkimagebasic
            milkkdtree
            milkinfo
            milkstatistic
            milkZernikePolyn
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
            $<INSTALL_INTERFACE:include>
            $<INSTALL_INTERFACE:include/${SRCNAME}>
    )
else()
    include("${CMAKE_CURRENT_SOURCE_DIR}/../cmake/coffee_standalone.cmake")
    add_library(${LIBNAME} SHARED ${SOURCEFILES})
    coffee_configure_standalone_target(${LIBNAME}
        MILK_LIBS
            milkimagegen
            milkfft
            milkimagebasic
            milkkdtree
            milkinfo
            milkstatistic
            milkZernikePolyn
            milkWFpropagate
            milkOpticsMaterials
    )
    install(TARGETS ${LIBNAME} DESTINATION lib)
endif()

install(FILES ${INCLUDEFILES} DESTINATION include/${SRCNAME})
install(PROGRAMS ${SCRIPTS} DESTINATION bin)

if(NOT TARGET coffee-all)
    add_custom_target(coffee-all
        DEPENDS
            coffeeOptSystProp
            coffeeAOsystSim
            coffeecoronagraphs
            coffeePIAACMCsimul
    )
endif()
