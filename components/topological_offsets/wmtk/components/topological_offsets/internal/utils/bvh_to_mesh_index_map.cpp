#include "bvh_to_mesh_index_map.hpp"

#include <wmtk/simplex/faces_single_dimension.hpp>
#include <wmtk/utils/Logger.hpp>
#include <wmtk/utils/TupleInspector.hpp>

namespace wmtk::components::internal::utils {

std::vector<int64_t> bvh_to_mesh_index_map(
    const SimpleBVH::BVH& bvh,
    const attribute::MeshAttributeHandle& position_handle)
{
    const Mesh& m = position_handle.mesh();
    const PrimitiveType pt_top = m.top_simplex_type();

    auto pos_acc = m.create_const_accessor<double>(position_handle);

    const auto faces = m.get_all(pt_top);

    std::vector<int64_t> index_map(faces.size(), -1);

    for (const Tuple& t : faces) {
        const auto vertices = m.orient_vertices(t);
        Eigen::VectorXd p_mid;
        p_mid.setZero(pos_acc.dimension());
        for (const Tuple& v : vertices) {
            p_mid += pos_acc.const_vector_attribute(v);
        }
        p_mid /= vertices.size();

        double sq_dist;
        SimpleBVH::VectorMax3d nearest_point;

        const int64_t nearest_facet = bvh.nearest_facet(p_mid, nearest_point, sq_dist);
        const int64_t fid = wmtk::utils::TupleInspector::global_cid(t);
        index_map[nearest_facet] = fid;
    }

    return index_map;
}

std::vector<int64_t> bvh_to_mesh_index_map(
    const SimpleBVH::BVH& bvh,
    const attribute::MeshAttributeHandle& position_handle,
    const PrimitiveType pt)
{
    const Mesh& m = position_handle.mesh();

    auto pos_acc = m.create_const_accessor<double>(position_handle);

    const auto faces = m.get_all(pt);

    std::vector<int64_t> index_map(faces.size(), -1);

    for (int64_t i = 0; i < faces.size(); ++i) {
        const Tuple& t = faces[i];

        const auto vertices = simplex::faces_single_dimension_tuples(
            m,
            simplex::Simplex(m, pt, t),
            PrimitiveType::Vertex);

        Eigen::VectorXd p_mid;
        p_mid.setZero(pos_acc.dimension());
        for (const Tuple& v : vertices) {
            p_mid += pos_acc.const_vector_attribute(v);
        }
        p_mid /= vertices.size();

        double sq_dist;
        SimpleBVH::VectorMax3d nearest_point;

        const int64_t nearest_facet = bvh.nearest_facet(p_mid, nearest_point, sq_dist);
        index_map[nearest_facet] = i;
    }

    bool missing_map = false;
    for (int64_t i = 0; i < index_map.size(); ++i) {
        if (index_map[i] == -1) {
            logger().error("No mapping for BVH index {}", i);
            missing_map = true;
        }
    }
    if (missing_map) {
        log_and_throw_error("Mapping from BVH to mesh is incomplete");
    }


    return index_map;
}

} // namespace wmtk::components::internal::utils