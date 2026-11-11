# =====================================================================
# coffeeAOsystSim - Adaptive Optics System Simulator
# =====================================================================

cmake_minimum_required(VERSION 3.16)

set(LIBNAME "coffeeAOsystSim")
set(SRCNAME "AOsystSim")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(SOURCEFILES
    AOsystSim.c
    AOsystSim_filter.c
    AOsystSim_extremeAO.c
    AOsystSim_pupil.c
    AOsystSim_wf_conf.c
    AOsystSim_wf.c
    AOsystSim_wfs_pyramid.c
    AOsystSim_pyrwfs.c
    AOsystSim_dm_shape.c
    AOsystSim_dm.c
    AOsystSim_lowfs_init.c
    AOsystSim_lowfs.c
    AOsystSim_fpwfs_imsimul.c
    AOsystSim_fpwfs_probes.c
    AOsystSim_fpwfs_solver.c
    AOsystSim_fpwfs_sens.c
    AOsystSim_cli.c
)

set(INCLUDEFILES
    AOsystSim.h
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
            milkinfo
            milkimagegen
            milkimagefilter
            milkstatistic
            milkfft
            coffeeOptSystProp
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
    if(TARGET coffeeOptSystProp)
        set(AOSIM_LOCAL_TARGETS LOCAL_TARGETS coffeeOptSystProp)
        set(AOSIM_EXTRA_MILK_LIBS "")
    else()
        set(AOSIM_LOCAL_TARGETS "")
        set(AOSIM_EXTRA_MILK_LIBS coffeeOptSystProp)
    endif()
    coffee_configure_standalone_target(${LIBNAME}
        MILK_LIBS
            milkinfo
            milkimagegen
            milkimagefilter
            milkstatistic
            milkfft
            milkWFpropagate
            milkOpticsMaterials
            ${AOSIM_EXTRA_MILK_LIBS}
        ${AOSIM_LOCAL_TARGETS}
    )
    install(TARGETS ${LIBNAME} DESTINATION lib)
endif()

install(FILES ${INCLUDEFILES} DESTINATION include/${SRCNAME})
