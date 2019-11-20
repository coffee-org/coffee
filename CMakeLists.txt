# run from _cmakebuild directory :
# cmake .. -DCMAKE_INSTALL_PREFIX=../_install
# cmake .. -DCMAKE_INSTALL_PREFIX=../_install -DUSE_MAGMA=ON
# cmake .. -DCMAKE_C_COMPILER=pgcc -DCMAKE_INSTALL_PREFIX=../_install -DUSE_MAGMA=ON
#
# make clean
# make
# sudo make install
#

cmake_minimum_required (VERSION 3.5)
project (coffee C)


# Version number
set ( VERSION_MAJOR 0 )
set ( VERSION_MINOR 1 )
set ( VERSION_PATCH 00 )
set ( VERSION_OPTION "-dev" )





# PERFORMANCE ORIENTED COMPILING OPTIONS

# Adds options to the compiler command line for targets in the current directory and below that are added after this command is invoked

add_compile_options(-std=gnu11)
add_compile_options(-Ofast)
add_compile_options(-march=native)

# runs the standard link-time optimizer
add_compile_options(-flto)

# A common flag is -pipe. This flag has no effect on the generated code, but it makes the compilation process faster. It tells the compiler to use pipes instead of temporary files during the different stages of compilation, which uses more memory. On systems with low memory, GCC might get killed. In those cases do not use this flag. 
add_compile_options(-pipe)


# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
message("Found OpenMP")
add_compile_options(-fopenmp)
endif()

message("COMPILE_OPTIONS:  ${COMPILE_OPTIONS}")




# options
option(USE_CUDA  "Use CUDA        library" OFF)
option(USE_MAGMA "Use MAGMA       library" OFF)
option(USE_MAGMA "Use QDWHpartial library" OFF)


# MAGMA (optional)
if(USE_MAGMA)
set(USE_CUDA ON)
endif(USE_MAGMA)


# QDWH (optional)
if(USE_QDWHPARTIAL)
set(USE_QDWHPARTIAL ON)
set(USE_CUDA ON)
set(USE_MAGMA ON)
endif(USE_QDWHPARTIAL)
set( QDWHPARTIAL_INCLUDE_PATH "/home/guyon/src/subaru/qdwhpartial-svd/magma" )
set( QDWHPARTIAL_LIBRARIES "psinv" )



if(USE_CUDA)
find_package( CUDA REQUIRED )
endif(USE_CUDA)

if(USE_MAGMA)
find_package(PkgConfig REQUIRED)
pkg_check_modules(MAGMA REQUIRED magma)
message(STATUS ${MAGMA_LIBRARY_DIRS})
link_directories( ${MAGMA_LIBRARY_DIRS} )
endif(USE_MAGMA)

# ncurses
SET(CURSES_USE_NCURSES TRUE)


message("====================================================")
message("PROJECT_SOURCE_DIR       = ${PROJECT_SOURCE_DIR}")

message("USE_CUDA                 = ${USE_CUDA}")

message("USE_MAGMA                = ${USE_MAGMA}")
if(USE_MAGMA)
message("   MAGMA_LIBRARIES       = ${MAGMA_LIBRARIES}")
endif(USE_MAGMA)

message("USE QDWHPARTIAL          = ${USE_QDWHPARTIAL}")
if(USE_QDWHPARTIAL)
message("   QDWHPARTIAL_LIBRARIES = ${QDWHPARTIAL_LIBRARIES}")
endif(USE_QDWHPARTIAL)

message("====================================================")




#                               
# recursively add libraries to be compiled
# each library compile as SHARED
#
add_subdirectory ( "src/ImageStreamIO" )
add_subdirectory ( "src/CommandLineInterface" )
add_subdirectory ( "src/00CORE" )
add_subdirectory ( "src/COREMOD_arith" )
add_subdirectory ( "src/COREMOD_iofits" )
add_subdirectory ( "src/COREMOD_memory" )
add_subdirectory ( "src/COREMOD_tools" )
add_subdirectory ( "src/info" )
add_subdirectory ( "src/fft" )
add_subdirectory ( "src/statistic" )
add_subdirectory ( "src/linopt_imtools" )
add_subdirectory ( "src/image_gen" )
add_subdirectory ( "src/image_filter" )
add_subdirectory ( "src/image_basic" )
add_subdirectory ( "src/ZernikePolyn" )
add_subdirectory ( "src/image_format" )
add_subdirectory ( "src/img_reduce" )
add_subdirectory ( "src/psf" )
add_subdirectory ( "src/cudacomp" )
add_subdirectory ( "src/kdtree" )
add_subdirectory ( "src/linARfilterPred" )

add_subdirectory ( "src/PIAACMCsimul" )
add_subdirectory ( "src/OptSystProp" )
add_subdirectory ( "src/OpticsMaterials" )
add_subdirectory ( "src/WFpropagate" )
add_subdirectory ( "src/coronagraphs" )
add_subdirectory ( "src/AOsystSim" )





# main
add_executable(coffee 
src/CLImain.c)

target_include_directories(coffee PUBLIC
${PROJECT_SOURCE_DIR}/src
${GSL_INCLUDE_DIRS}
${FFTW_INCLUDE_DIRS}
${FFTWF_INCLUDE_DIRS}
${CURSES_INCLUDE_DIR}
)

# adds the options to all targets within the directory and its sub-directories
target_compile_options(coffee PUBLIC 
${FFTW_CFLAGS_OTHER}
${FFTWF_CFLAGS_OTHER}
-fopenmp
)

set(BUILD_FLAGS "-DPACKAGE_NAME=\\\"coffee\\\" -DPACKAGE_VERSION=\\\"1.0.00\\\" -DCONFIGDIR=\\\"${PROJECT_SOURCE_DIR}/config\\\" -DSOURCEDIR=\\\"${PROJECT_SOURCE_DIR}/src\\\" -DABSSRCTOPDIR=\\\"${PROJECT_SOURCE_DIR}\\\" -DPACKAGE_BUGREPORT=\\\"https://github.com/coffee-org/coffee/issues\\\"")


target_link_libraries (coffee PUBLIC
CLIcore 
00CORE 
COREMOD_arith 
COREMOD_memory 
ImageStreamIO 
COREMOD_tools 
info 
fft
linopt_imtools 
image_basic
COREMOD_iofits 
ZernikePolyn
image_filter
image_gen
statistic
image_format 
img_reduce 
psf 
cudacomp 
kdtree 
linARfilterPred
PIAACMCsimul
OptSystProp
OpticsMaterials
WFpropagate
coronagraphs
AOsystSim
m
readline
ncurses
cfitsio
dl
rt
${GSL_LIBRARIES} 
${FFTW_LIBRARIES}
${FFTWF_LIBRARIES}
${CMAKE_THREAD_LIBS_INIT}
-fopenmp
)




if(USE_MAGMA)
target_link_libraries (coffee PUBLIC ${MAGMA_LIBRARIES})
set(BUILD_FLAGS "${BUILD_FLAGS} -DHAVE_CUDA -DHAVE_MAGMA" )
endif(USE_MAGMA)

if(USE_QDWHPARTIAL)
target_link_libraries (coffee PUBLIC ${QDWHPARTIAL_LIBRARIES})
target_include_directories(coffee PUBLIC ${QDWHPARTIAL_INCLUDE_PATH} )
set(BUILD_FLAGS "${BUILD_FLAGS} -DHAVE_QDWHpartial" )
endif(USE_QDWHPARTIAL)


set_target_properties(coffee PROPERTIES COMPILE_FLAGS "${BUILD_FLAGS}" )

message("====================================================")
get_target_property(coffeecompopts coffee COMPILE_OPTIONS)
message("BUILD_FLAGS    :  ${BUILD_FLAGS}")
message("COMPILE_OPTIONS:  ${coffeecompopts}")
message("CMAKE_EXE_LINKER_FLAGS  :  ${CMAKE_EXE_LINKER_FLAGS}")
message("====================================================")

install(TARGETS coffee DESTINATION bin)
