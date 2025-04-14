include(wmtk_add_application)
wmtk_add_application(topological_offsets_app
    topological_offsets_main.cpp
    topological_offsets_spec.hpp
    utils/igl_input.hpp
    utils/igl_input.cpp
    )


# topological_offsets requires the input component and the topological_offsets component
target_link_libraries(topological_offsets_app PRIVATE
wmtk::input
wmtk::topological_offsets
wmtk::output 
wmtk::multimesh_from_tag
nlohmann_json::nlohmann_json
)