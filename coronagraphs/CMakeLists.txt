# library name
set(LIBNAME "coffeecoronagraphs")
set(SRCNAME "coronagraphs")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(SOURCEFILES
    coronagraphs.c
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
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/coremods
        ${CMAKE_SOURCE_DIR}/plugins/milk-extra-src
        ${CMAKE_SOURCE_DIR}/plugins/milkatmturb/src
        $<TARGET_PROPERTY:CLIcore,INTERFACE_INCLUDE_DIRECTORIES>
        /usr/local/milk/include
    )

    set(OPT_WF_LIBS "")
    if(TARGET milkWFpropagate)
        list(APPEND OPT_WF_LIBS milkWFpropagate)
    elseif(TARGET milkatmturb)
        list(APPEND OPT_WF_LIBS milkatmturb)
    else()
        list(APPEND OPT_WF_LIBS milkWFpropagate)
    endif()

    target_link_libraries(${LIBNAME} PUBLIC
        CLIcore
        ImageStreamIO
        milkCOREMODmemory
        milkCOREMODarith
        milkCOREMODiofits
        milkCOREMODtools
        milkimagegen
        milkfft
        milkimagebasic
        milkinfo
        milkstatistic
        milkZernikePolyn
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
            /usr/local/milk/include
            /usr/local/include
        )

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
    endif()

    target_include_directories(${LIBNAME} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${MILK_STANDALONE_INCLUDES}
    )

    set(MILK_LIB_NAMES
        CLIcore
        ImageStreamIO
        milkCOREMODmemory
        milkCOREMODarith
        milkCOREMODiofits
        milkCOREMODtools
        milkimagegen
        milkfft
        milkimagebasic
        milkinfo
        milkstatistic
        milkZernikePolyn
        milkWFpropagate
    )

    set(MILK_SEARCH_HINTS
        ${MILK_ROOT}/lib
        $ENV{MILK_ROOT}/lib
        /usr/local/milk/lib
        /usr/local/lib
    )
    if(MILK_SOURCE_FOUND)
        list(APPEND MILK_SEARCH_HINTS
            ${MILK_SOURCE_FOUND}/_build/lib
            ${MILK_SOURCE_FOUND}/_build/src/cli/CLIcore
            ${MILK_SOURCE_FOUND}/_build/src/engine/ImageStreamIO
            ${MILK_SOURCE_FOUND}/_build/plugins/milkatmturb
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/fft
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/image_gen
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/image_basic
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/info
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/statistic
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/ZernikePolyn
        )
    endif()

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

if(NOT TARGET coffee-all)
    add_custom_target(coffee-all
        DEPENDS
            coffeeOptSystProp
            coffeeAOsystSim
            coffeecoronagraphs
            coffeePIAACMCsimul
    )
endif()

