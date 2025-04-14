#pragma once

#include <SimpleBVH/BVH.hpp>
#include <wmtk/Mesh.hpp>

#include "OffsetBvh.hpp"

namespace wmtk::components::internal::utils {
/**
 * @brief Check if tet is fully covered by the offset.
 *
 * The method samples the tet and is therefore not exact. The current sampling just uses 5 points,
 * the corners and the center. The sampling is subject to change.
 */
bool is_tet_in_offset(
    const attribute::MeshAttributeHandle& pos_handle,
    const SimpleBVH::BVH& bvh,
    const double offset_distance,
    const simplex::Simplex& s);

bool is_tet_in_offset(
    const attribute::MeshAttributeHandle& pos_handle,
    const OffsetBvh& bvh,
    const simplex::Simplex& s);

/**
 * @brief A conservative check for a tet being inside the offset.
 *
 * The method covers space with voxels. Each voxel is approximated by a sphere. First, every voxel
 * that does not intersect the tet is removed from the queue. If all remaining voxels are within the
 * offset, the entire tet must be. If one voxel is outside, the entire tet is considered outside.
 *
 * If a voxel is intersected by the offset, it is subdivided into 8 smaller voxels and those are
 * added to the queue. This hierarchical approach is repeated until voxels are smaller than a given
 * threshold. If any voxel reaches that threshold, we conservatively assume that the tet is outside.
 */
bool is_tet_in_offset_conservative_sampling(
    const Eigen::MatrixXd& bbox,
    const attribute::MeshAttributeHandle& pos_handle,
    const OffsetBvh& bvh,
    const simplex::Simplex& s,
    const double threshold);

/**
 * @brief An aggressive check for a tet bein inside the offset.
 *
 * Just like the conservative but returns true if any of the samples is inside.
 */
bool is_tet_in_offset_aggressive_sampling(
    const Eigen::MatrixXd& bbox,
    const attribute::MeshAttributeHandle& pos_handle,
    const OffsetBvh& bvh,
    const simplex::Simplex& s,
    const double threshold);

/**
 * @brief A conservative estimate if a tet is inside the offset.
 *
 * The tet is approximated by a sphere.
 */
bool is_tet_in_offset_conservative(
    const attribute::MeshAttributeHandle& pos_handle,
    const SimpleBVH::BVH& bvh,
    const double offset_distance,
    const simplex::Simplex& s);

} // namespace wmtk::components::internal::utils
