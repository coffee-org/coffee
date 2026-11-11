# =====================================================================
# coffeeAOsystSim - Adaptive Optics System Simulator
# =====================================================================

set(LIBNAME "coffeeAOsystSim")
set(SRCNAME "AOsystSim")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(SOURCEFILES
    AOsystSim.c
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
    coffee_configure_standalone_target(${LIBNAME}
        MILK_LIBS
            milkinfo
            milkimagegen
            milkimagefilter
            milkstatistic
            milkfft
            milkWFpropagate
            milkOpticsMaterials
        LOCAL_TARGETS
            coffeeOptSystProp
    )
    install(TARGETS ${LIBNAME} DESTINATION lib)
endif()

install(FILES ${INCLUDEFILES} DESTINATION include/${SRCNAME})
