#pragma once

#include <wmtk/Mesh.hpp>

namespace wmtk::utils {

/**
 * @brief Compute the bounding box of a mesh.
 *
 * @param position_handle The position handle, which also contains a pointer to the mesh.
 * @return A matrix where the first row contains all minimal coordinates and the second row contains
 * all maximal coordinates.
 */
Eigen::MatrixXd bbox_from_mesh(const attribute::MeshAttributeHandle& position_handle);

/**
 * @brief Compute the bounding box of a tagged part of the mesh.
 *
 * @param position_handle The position handle, which also contains a pointer to the mesh.
 * @param position_handle The tag handle.
 * @return A matrix where the first row contains all minimal coordinates and the second row contains
 * all maximal coordinates.
 */
Eigen::MatrixXd bbox_from_mesh(
    const attribute::MeshAttributeHandle& position_handle,
    const attribute::MeshAttributeHandle& tag_handle,
    const int64_t tag_value);

/**
 * @brief Compute the diagonal of the bounding box of a mesh.
 * @param position_handle The position handle, which also contains a pointer to the mesh.
 * @return The diagonal of the bounding box.
 */
double bbox_diagonal_from_mesh(const attribute::MeshAttributeHandle& position_handle);

double bbox_diagonal_from_mesh(
    const attribute::MeshAttributeHandle& position_handle,
    const attribute::MeshAttributeHandle& tag_handle,
    const int64_t tag_value);

} // namespace wmtk::utils