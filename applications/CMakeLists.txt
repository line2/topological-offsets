list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/applications/cmake")


# clear out the parent targets
set(WMTK_APPLICATION_HEADER_TARGETS )

include(wmtk_add_application)
include(wmtk_generate_test_config)

option(WMTK_ENABLE_APPLICATION_topological_offsets "Runs topological_offsets application" ON)
if(WMTK_ENABLE_APPLICATION_topological_offsets)
add_subdirectory(topological_offsets)
add_subdirectory(image_simulation)
endif()