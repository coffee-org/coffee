project(lib_AOsystSim_project)

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${PROJECT_SOURCE_DIR}/..")
add_library(AOsystSim SHARED AOsystSim.c)

install(TARGETS AOsystSim DESTINATION lib)
install(FILES AOsystSim.h DESTINATION include)

