#pragma once

#include <igl/Hit.h>
#include <wmtk/Mesh.hpp>

namespace wmtk::components::internal::utils {

/**
 * @brief A slow wrapper of the even slower ray_mesh_intersect from IGL.
 */
bool ray_mesh_intersect(
    const attribute::MeshAttributeHandle& position_handle,
    const Eigen::Vector3d& source,
    const Eigen::Vector3d& dir,
    igl::Hit& hit);

} // namespace wmtk::components::internal::utils