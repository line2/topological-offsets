
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/components/cmake/recipes")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/components/cmake")

include(add_component)

# wmtk::component_utils
add_subdirectory(utils)

# list of all components
add_subdirectory("input")
add_subdirectory("multimesh")
add_subdirectory("output/wmtk/components/output")
add_subdirectory("multimesh_from_tag/wmtk/components/multimesh_from_tag")
add_subdirectory("marching/wmtk/components/marching")
add_subdirectory("regular_space/wmtk/components/regular_space")
add_subdirectory("topological_offsets/wmtk/components/topological_offsets")






