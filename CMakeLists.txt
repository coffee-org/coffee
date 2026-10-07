# =====================================================================
# coffee - Coronagraph Optimization For Fast Exoplanets Exploration
# =====================================================================

cmake_minimum_required(VERSION 3.16)

# In-tree build detection:
# If CLIcore target already exists, we are building inside the milk framework.
# Milk's plugin discovery scanner automatically discovers and adds each module
# subdirectory (OptSystProp, AOsystSim, PIAACMCsimul, coronagraphs) at depth 2.
# Return immediately to avoid duplicate subdirectory additions.
if(TARGET CLIcore)
    return()
endif()

project(coffee C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Find OpenMP (required for multi-threaded simulations)
find_package(OpenMP REQUIRED)

# Find GSL (GNU Scientific Library)
find_package(GSL REQUIRED)

# Find FFTW3
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED fftw3)
pkg_check_modules(FFTWF REQUIRED fftw3f)

# Discover MILK framework source and installation paths
set(MILK_SOURCE_CANDIDATES
    ${MILK_SOURCE_DIR}
    $ENV{MILK_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../milk-framework-dev
    $ENV{MILK_ROOT}
    ${CMAKE_CURRENT_SOURCE_DIR}/../milk
    /usr/local/milk
)

set(MILK_SOURCE_FOUND "")
foreach(cand ${MILK_SOURCE_CANDIDATES})
    if(cand AND EXISTS "${cand}/src/coremods")
        get_filename_component(MILK_SOURCE_FOUND "${cand}" ABSOLUTE)
        message(STATUS "Found MILK source tree: ${MILK_SOURCE_FOUND}")
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

set(MILK_STANDALONE_INCLUDES
    ${CMAKE_CURRENT_SOURCE_DIR}
)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../AtmosphericTurbulenceSimul/src")
    list(APPEND MILK_STANDALONE_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/../AtmosphericTurbulenceSimul/src")
endif()

if(MILK_SOURCE_FOUND)
    list(APPEND MILK_STANDALONE_INCLUDES
        ${MILK_SOURCE_FOUND}/src
        ${MILK_SOURCE_FOUND}/src/coremods
        ${MILK_SOURCE_FOUND}/plugins/milk-extra-src
        ${MILK_SOURCE_FOUND}/plugins/milkatmturb/src
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

list(APPEND MILK_STANDALONE_INCLUDES
    /usr/local/milk/include
    /usr/local/include
)

# Add subdirectories
add_subdirectory(OptSystProp)
add_subdirectory(AOsystSim)
add_subdirectory(coronagraphs)
add_subdirectory(PIAACMCsimul)

# Top-level meta target
if(NOT TARGET coffee-all)
    add_custom_target(coffee-all ALL
        DEPENDS
            coffeeOptSystProp
            coffeeAOsystSim
            coffeecoronagraphs
            coffeePIAACMCsimul
    )
endif()

