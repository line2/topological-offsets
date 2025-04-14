#pragma once

#include <SimpleBVH/BVH.hpp>
#include <memory>
#include <utility>
#include <wmtk/Mesh.hpp>

namespace wmtk::components::internal::utils {

class OffsetBvh
{
public:
    OffsetBvh(const attribute::MeshAttributeHandle& position_handle);
    OffsetBvh(const attribute::MeshAttributeHandle& position_handle, const PrimitiveType pt);

    void update_face_id_map();

    double sq_dist(const SimpleBVH::VectorMax3d& p, SimpleBVH::VectorMax3d& nearest_point) const;

    double sq_dist(const SimpleBVH::VectorMax3d& p) const;

    double dist(const SimpleBVH::VectorMax3d& p) const;

    /**
     * @brief Returns the squared distance and the Tuple of the nearest facet.
     */
    std::tuple<double, Tuple> sq_dist_and_tuple(
        const SimpleBVH::VectorMax3d& p,
        SimpleBVH::VectorMax3d& nearest_point) const;

    std::tuple<double, Tuple> sq_dist_and_tuple(const SimpleBVH::VectorMax3d& p) const;

    std::tuple<double, double> sq_dist_and_offset_distance(
        const SimpleBVH::VectorMax3d& p,
        SimpleBVH::VectorMax3d& nearest_point) const;

    std::tuple<double, double> sq_dist_and_offset_distance(const SimpleBVH::VectorMax3d& p) const;

    std::tuple<double, double, SimpleBVH::VectorMax3d> sq_dist_offset_distance_nearest_point(
        const SimpleBVH::VectorMax3d& p) const;

    attribute::MeshAttributeHandle f_offset_distance_handle() const;
    attribute::MeshAttributeHandle v_offset_distance_handle() const;
    void register_offset_distance_attribute(const double offset_distance);

    void update_offset_distance_handle();
    void update_position_handle();

    const SimpleBVH::BVH& bvh() const;

    std::shared_ptr<SimpleBVH::BVH> bvh_ptr() const;

    Mesh& mesh();
    const Mesh& mesh() const;

private:
    PrimitiveType m_pt;
    std::shared_ptr<SimpleBVH::BVH> m_bvh;

    attribute::MeshAttributeHandle m_position_handle;
    attribute::MeshAttributeHandle m_f_offset_distance_handle;
    attribute::MeshAttributeHandle m_v_offset_distance_handle;

    std::vector<int64_t> m_bvh_to_mesh_index;

    const std::string m_f_offset_distance_handle_name = "f_offset_bvh_offset_distance";
    const std::string m_v_offset_distance_handle_name = "v_offset_bvh_offset_distance";
};

} // namespace wmtk::components::internal::utils