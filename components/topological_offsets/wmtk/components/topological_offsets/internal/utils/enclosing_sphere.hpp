#pragma once

#include <wmtk/Mesh.hpp>

namespace wmtk::components::internal::utils {

/**
 * @brief Compute a sphere that is enclosing the input simplex.
 *
 * The sphere is computed by taking the center and finding the vertex that is the furthest away from
 * it.
 *
 * @param pos_handle The position handle of the mesh.
 * @param s The simplex for which the enclosing sphere should be computed.
 *
 * @returns The center point and the sphere radius.
 */
std::tuple<Eigen::VectorXd, double> enclosing_sphere(
    const attribute::MeshAttributeHandle& pos_handle,
    const simplex::Simplex& s);

} // namespace wmtk::components::internal::utils
