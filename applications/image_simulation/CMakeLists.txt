include(wmtk_add_application)
wmtk_add_application(image_simulation_app
    image_simulation_main.cpp
    image_simulation_spec.hpp
    utils/igl_input.hpp
    utils/igl_input.cpp
    utils/write_msh.hpp
    )


# image_simulation requires the input component and the image_simulation component
target_link_libraries(image_simulation_app PRIVATE
wmtk::input
wmtk::topological_offsets
wmtk::output 
wmtk::multimesh_from_tag
nlohmann_json::nlohmann_json
)