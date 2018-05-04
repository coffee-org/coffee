project(lib_coronagraphs_project)

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${PROJECT_SOURCE_DIR}/..")
add_library(coronagraphs SHARED coronagraphs.c)

install(TARGETS coronagraphs DESTINATION lib)
install(FILES coronagraphs.h DESTINATION include)

