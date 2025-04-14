#pragma once

#include <wmtk/attribute/MeshAttributeHandle.hpp>

namespace wmtk::components::utils {

/*
 * @brief Write mesh in Paraview file format.
 *
 * This method is supposed to help in debugging.
 */
void write_mesh(
    const attribute::MeshAttributeHandle& pos_handle,
    const std::string& name,
    const bool intermediate_output = true);

void write_mesh(
    const attribute::MeshAttributeHandle& pos_handle,
    const std::string& name,
    const PrimitiveType& pt,
    const bool intermediate_output = true);

} // namespace wmtk::components::utils