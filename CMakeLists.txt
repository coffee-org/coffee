# =====================================================================
# coffee - Coronagraph Optimization For Fast Exoplanets Exploration
# =====================================================================

cmake_minimum_required(VERSION 3.16)

# In-tree build detection:
# When building inside the milk framework, milk's plugin discovery scanner
# automatically discovers and adds each module subdirectory at depth 2.
if(COMMAND milk_add_runnable_module OR TARGET CLIcore)
    return()
endif()

project(coffee C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED fftw3)
pkg_check_modules(FFTWF REQUIRED fftw3f)

# Discover MILK framework source and installation paths
set(MILK_SOURCE_CANDIDATES
    ${MILK_SOURCE_DIR}
    $ENV{MILK_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src-frameworkdev/milk
    ${CMAKE_CURRENT_SOURCE_DIR}/../milk
    ${CMAKE_CURRENT_SOURCE_DIR}/../milk-framework-dev
    $ENV{MILK_ROOT}
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

set(MILK_ROOT_CANDIDATES
    ${MILK_ROOT}
    $ENV{MILK_ROOT}
    $ENV{MILK_INSTALLDIR}
)
if(MILK_SOURCE_FOUND)
    file(GLOB _milk_ver_installs "${MILK_SOURCE_FOUND}/_install/milk-*")
    list(APPEND MILK_ROOT_CANDIDATES ${_milk_ver_installs} "${MILK_SOURCE_FOUND}/_install")
endif()
if(NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    list(APPEND MILK_ROOT_CANDIDATES "${CMAKE_INSTALL_PREFIX}")
endif()
list(APPEND MILK_ROOT_CANDIDATES /usr/local/milk /usr/local)

set(MILK_ROOT_FOUND "")
foreach(rcand ${MILK_ROOT_CANDIDATES})
    if(rcand AND EXISTS "${rcand}/include/libcomputeunit/milk_runnable.h")
        get_filename_component(MILK_ROOT_FOUND "${rcand}" ABSOLUTE)
        message(STATUS "Found MILK install root: ${MILK_ROOT_FOUND}")
        break()
    endif()
endforeach()

if(MILK_ROOT_FOUND AND NOT MILK_ROOT)
    set(MILK_ROOT "${MILK_ROOT_FOUND}")
endif()

find_path(MILK_INCLUDE_DIR
    NAMES libcomputeunit/milk_runnable.h
    HINTS
        ${MILK_ROOT}/include
        $ENV{MILK_ROOT}/include
        $ENV{MILK_INSTALLDIR}/include
        ${CMAKE_INSTALL_PREFIX}/include
    NO_DEFAULT_PATH
)
if(NOT MILK_INCLUDE_DIR)
    find_path(MILK_INCLUDE_DIR
        NAMES libcomputeunit/milk_runnable.h CommandLineInterface/CLIcore.h
        HINTS /usr/local/milk/include /usr/local/include
    )
endif()

set(MILK_STANDALONE_INCLUDES "")

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../OpticsMaterials")
    list(APPEND MILK_STANDALONE_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/../OpticsMaterials")
endif()
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../WFpropagate")
    list(APPEND MILK_STANDALONE_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/../WFpropagate")
endif()

if(MILK_INCLUDE_DIR)
    list(APPEND MILK_STANDALONE_INCLUDES
        ${MILK_INCLUDE_DIR}
        ${MILK_INCLUDE_DIR}/CommandLineInterface
        ${MILK_INCLUDE_DIR}/CommandLineInterface/CLIcore
        ${MILK_INCLUDE_DIR}/libmilkcommon
        ${MILK_INCLUDE_DIR}/libmilkdata
        ${MILK_INCLUDE_DIR}/libfps
        ${MILK_INCLUDE_DIR}/libprocessinfo
        ${MILK_INCLUDE_DIR}/libcomputeunit
        ${MILK_INCLUDE_DIR}/ImageStreamIO
    )
endif()

if(MILK_SOURCE_FOUND)
    list(APPEND MILK_STANDALONE_INCLUDES
        ${MILK_SOURCE_FOUND}/src
        ${MILK_SOURCE_FOUND}/src/coremods
        ${MILK_SOURCE_FOUND}/plugins/milk-extra-src
        ${MILK_SOURCE_FOUND}/plugins/OpticsMaterials
        ${MILK_SOURCE_FOUND}/plugins/WFpropagate
        ${MILK_SOURCE_FOUND}/src/cli
        ${MILK_SOURCE_FOUND}/src/cli/CLIcore
        ${MILK_SOURCE_FOUND}/src/cli/libmilkscript
        ${MILK_SOURCE_FOUND}/src/engine
        ${MILK_SOURCE_FOUND}/src/engine/libfps
        ${MILK_SOURCE_FOUND}/src/engine/libmilkcommon
        ${MILK_SOURCE_FOUND}/src/engine/libmilkdata
        ${MILK_SOURCE_FOUND}/src/engine/libprocessinfo
        ${MILK_SOURCE_FOUND}/src/engine/libcomputeunit
        ${MILK_SOURCE_FOUND}/src/engine/ImageStreamIO
    )
    if(EXISTS "${MILK_SOURCE_FOUND}/_build/config")
        list(APPEND MILK_STANDALONE_INCLUDES "${MILK_SOURCE_FOUND}/_build/config")
    endif()
endif()

set(MILK_SEARCH_HINTS "")
if(MILK_ROOT)
    list(APPEND MILK_SEARCH_HINTS "${MILK_ROOT}/lib")
endif()
if(DEFINED ENV{MILK_INSTALLDIR})
    list(APPEND MILK_SEARCH_HINTS "$ENV{MILK_INSTALLDIR}/lib")
endif()
list(APPEND MILK_SEARCH_HINTS
    ${CMAKE_INSTALL_PREFIX}/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/../OpticsMaterials/_build
    ${CMAKE_CURRENT_SOURCE_DIR}/../WFpropagate/_build
)
if(MILK_SOURCE_FOUND)
    list(APPEND MILK_SEARCH_HINTS
        ${MILK_SOURCE_FOUND}/_build/lib
        ${MILK_SOURCE_FOUND}/_build/src/engine/libcomputeunit
        ${MILK_SOURCE_FOUND}/_build/src/engine/libfps
        ${MILK_SOURCE_FOUND}/_build/src/engine/libprocessinfo
        ${MILK_SOURCE_FOUND}/_build/src/engine/libmilkdata
        ${MILK_SOURCE_FOUND}/_build/src/engine/ImageStreamIO
        ${MILK_SOURCE_FOUND}/_build/src/coremods/COREMOD_memory
        ${MILK_SOURCE_FOUND}/_build/src/coremods/COREMOD_arith
        ${MILK_SOURCE_FOUND}/_build/src/coremods/COREMOD_iofits
        ${MILK_SOURCE_FOUND}/_build/src/coremods/COREMOD_tools
        ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/fft
        ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/image_gen
        ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/image_basic
        ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/image_filter
        ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/info
        ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/statistic
        ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/ZernikePolyn
        ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/linopt_imtools
    )
endif()
list(APPEND MILK_SEARCH_HINTS /usr/local/milk/lib /usr/local/lib)

# Add subdirectories in topological dependency order
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

# Refresh milk command catalog upon out-of-tree installation into a milk prefix
install(CODE "
    set(_milk_prefix \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}\")
    if(EXISTS \"\${_milk_prefix}/bin/milk-cli\")
        execute_process(COMMAND \${CMAKE_COMMAND} -E env
                                MILK_INSTALLDIR=\${_milk_prefix}
                                LD_LIBRARY_PATH=\${_milk_prefix}/lib:\$ENV{LD_LIBRARY_PATH}
                                \${_milk_prefix}/bin/milk-cli --catalog
                        RESULT_VARIABLE _milk_rc)
    endif()
")
