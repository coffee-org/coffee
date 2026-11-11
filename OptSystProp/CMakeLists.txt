# =====================================================================
# coffeeOptSystProp - Optical System Fresnel Propagation Engine
# =====================================================================

set(LIBNAME "coffeeOptSystProp")
set(SRCNAME "OptSystProp")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(SOURCEFILES
    OptSystProp.c
    OptSystProp_propagateCube.c
    OptSystProp_run.c
)

set(INCLUDEFILES
    OptSystProp.h
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
            milkfft
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
            milkWFpropagate
            milkOpticsMaterials
            milkfft
    )
    install(TARGETS ${LIBNAME} DESTINATION lib)
endif()

install(FILES ${INCLUDEFILES} DESTINATION include/${SRCNAME})
