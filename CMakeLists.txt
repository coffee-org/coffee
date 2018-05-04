project(lib_OptSystProp_project)

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${PROJECT_SOURCE_DIR}/..")


add_library(OptSystProp SHARED OptSystProp.h OptSystProp.c  OptSystProp_propagateCube.c  OptSystProp_run.c )


install(TARGETS OptSystProp DESTINATION lib)
install(FILES OptSystProp.h DESTINATION include)

