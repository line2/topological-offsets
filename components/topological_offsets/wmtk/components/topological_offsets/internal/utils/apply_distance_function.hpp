#pragma once

#include <wmtk/Mesh.hpp>

namespace wmtk::components::internal::utils {

/**
 * @brief Evaluate the tinyexpr expression on every vertex and write the value into the (vertex)
 * distance handle.
 */
void apply_distance_function(
    attribute::MeshAttributeHandle& pos_handle,
    attribute::MeshAttributeHandle& distance_handle,
    double distance,
    const std::string& expression);

} // namespace wmtk::components::internal::utils