# =====================================================================
# coffeecoronagraphs - Coronagraph Design & Propagation Algorithms
# =====================================================================

set(LIBNAME "coffeecoronagraphs")
set(SRCNAME "coronagraphs")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(SOURCEFILES
    coronagraphs.c
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
