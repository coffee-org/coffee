project(lib_coronagraphs_project)

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${PROJECT_SOURCE_DIR}/..")
add_library(coronagraphs SHARED coronagraphs.c)

set_target_properties(coronagraphs PROPERTIES COMPILE_FLAGS "-DCORONAGRAPHSDATADIR=\\\"${PROJECT_SOURCE_DIR}/src/coronagraphs/datadir\\\"" )


install(TARGETS coronagraphs DESTINATION lib)
install(FILES coronagraphs.h DESTINATION include)

