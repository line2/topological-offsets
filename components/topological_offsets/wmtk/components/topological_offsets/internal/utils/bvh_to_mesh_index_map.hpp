#pragma once

#include <SimpleBVH/BVH.hpp>
#include <vector>
#include <wmtk/Mesh.hpp>

namespace wmtk::components::internal::utils {

/**
 * @brief Create a map from the bvh nearest facet index to the mesh face index.
 *
 * The bvh must be created for the same mesh that is given by the position handle.
 */
std::vector<int64_t> bvh_to_mesh_index_map(
    const SimpleBVH::BVH& bvh,
    const attribute::MeshAttributeHandle& position_handle);

std::vector<int64_t> bvh_to_mesh_index_map(
    const SimpleBVH::BVH& bvh,
    const attribute::MeshAttributeHandle& position_handle,
    const PrimitiveType pt);

} // namespace wmtk::components::internal::utils