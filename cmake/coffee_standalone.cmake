# =====================================================================
# coffee_standalone.cmake - Shared out-of-tree MILK discovery & linkage
# =====================================================================

function(coffee_configure_standalone_target TARGET_NAME)
    cmake_parse_arguments(ARG "" "" "EXTRA_INCLUDES;MILK_LIBS;LOCAL_TARGETS" ${ARGN})

    target_compile_definitions(${TARGET_NAME} PRIVATE MILK_NO_CLI)
    target_compile_options(${TARGET_NAME} PRIVATE
        -Wall -Wextra -Wstrict-prototypes -Wno-unused-const-variable
    )

    if(NOT DEFINED MILK_STANDALONE_INCLUDES)
        set(_coffee_root "${CMAKE_CURRENT_LIST_DIR}/..")
        set(_src_cands
            ${MILK_SOURCE_DIR}
            $ENV{MILK_SOURCE_DIR}
            "${_coffee_root}/../../src-frameworkdev/milk"
            "${_coffee_root}/../milk"
            "${_coffee_root}/../milk-framework-dev"
            $ENV{MILK_ROOT}
            /usr/local/milk
        )
        set(MILK_SOURCE_FOUND "")
        foreach(cand ${_src_cands})
            if(cand AND EXISTS "${cand}/src/coremods")
                get_filename_component(MILK_SOURCE_FOUND "${cand}" ABSOLUTE)
                break()
            endif()
        endforeach()

        set(_root_cands ${MILK_ROOT} $ENV{MILK_ROOT} $ENV{MILK_INSTALLDIR})
        if(MILK_SOURCE_FOUND)
            file(GLOB _ver_inst "${MILK_SOURCE_FOUND}/_install/milk-*")
            list(APPEND _root_cands ${_ver_inst} "${MILK_SOURCE_FOUND}/_install")
        endif()
        list(APPEND _root_cands /usr/local/milk /usr/local)
        foreach(rcand ${_root_cands})
            if(rcand AND EXISTS "${rcand}/include/libcomputeunit/milk_runnable.h")
                get_filename_component(MILK_ROOT "${rcand}" ABSOLUTE)
                break()
            endif()
        endforeach()

        find_path(MILK_INCLUDE_DIR
            NAMES libcomputeunit/milk_runnable.h CommandLineInterface/CLIcore.h
            HINTS ${MILK_ROOT}/include /usr/local/milk/include /usr/local/include
        )

        set(MILK_STANDALONE_INCLUDES "")
        if(EXISTS "${_coffee_root}/../OpticsMaterials")
            list(APPEND MILK_STANDALONE_INCLUDES "${_coffee_root}/../OpticsMaterials")
        endif()
        if(EXISTS "${_coffee_root}/../WFpropagate")
            list(APPEND MILK_STANDALONE_INCLUDES "${_coffee_root}/../WFpropagate")
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
                ${MILK_SOURCE_FOUND}/src/cli
                ${MILK_SOURCE_FOUND}/src/cli/CLIcore
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

        set(MILK_SEARCH_HINTS
            ${MILK_ROOT}/lib
            ${CMAKE_INSTALL_PREFIX}/lib
            "${_coffee_root}/../OpticsMaterials/_build"
            "${_coffee_root}/../WFpropagate/_build"
            /usr/local/milk/lib
            /usr/local/lib
        )
    endif()

    target_include_directories(${TARGET_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
    )
    if(ARG_EXTRA_INCLUDES)
        target_include_directories(${TARGET_NAME} PRIVATE ${ARG_EXTRA_INCLUDES})
    endif()
    target_include_directories(${TARGET_NAME} SYSTEM PRIVATE ${MILK_STANDALONE_INCLUDES})

    set(_base_libs
        milkcomputeunit
        milkfps
        milkprocessinfo
        milkdata
        ImageStreamIO
        milkCOREMODmemory
        milkCOREMODarith
        milkCOREMODiofits
        milkCOREMODtools
        ${ARG_MILK_LIBS}
    )

    set(_resolved_libs "")
    foreach(lib ${_base_libs})
        find_library(_COFFEE_LIB_${lib} NAMES ${lib} HINTS ${MILK_SEARCH_HINTS} NO_DEFAULT_PATH)
        if(NOT _COFFEE_LIB_${lib})
            find_library(_COFFEE_LIB_${lib} NAMES ${lib} HINTS ${MILK_SEARCH_HINTS})
        endif()
        if(_COFFEE_LIB_${lib})
            list(APPEND _resolved_libs ${_COFFEE_LIB_${lib}})
        else()
            list(APPEND _resolved_libs ${lib})
        endif()
    endforeach()

    set(_rpath "$ORIGIN")
    if(MILK_ROOT)
        list(APPEND _rpath "${MILK_ROOT}/lib")
    endif()
    set_target_properties(${TARGET_NAME} PROPERTIES
        LINKER_LANGUAGE C
        INTERPROCEDURAL_OPTIMIZATION OFF
        INSTALL_RPATH "${_rpath}"
    )

    target_link_libraries(${TARGET_NAME} PUBLIC
        ${_resolved_libs}
        ${ARG_LOCAL_TARGETS}
        OpenMP::OpenMP_C
        GSL::gsl
        GSL::gslcblas
        ${FFTW_LIBRARIES}
        ${FFTWF_LIBRARIES}
        m
    )
endfunction()
