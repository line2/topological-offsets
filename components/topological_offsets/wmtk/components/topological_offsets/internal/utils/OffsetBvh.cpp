#include "OffsetBvh.hpp"

#include "bvh_from_mesh.hpp"
#include "bvh_to_mesh_index_map.hpp"

namespace wmtk::components::internal::utils {

OffsetBvh::OffsetBvh(const attribute::MeshAttributeHandle& position_handle)
    : m_position_handle(position_handle)
    , m_bvh(bvh_from_mesh(position_handle))
    , m_pt(position_handle.mesh().top_simplex_type())
{
    update_face_id_map();
}

OffsetBvh::OffsetBvh(const attribute::MeshAttributeHandle& position_handle, const PrimitiveType pt)
    : m_position_handle(position_handle)
    , m_bvh(bvh_from_mesh(position_handle, pt))
    , m_pt(pt)
{
    update_face_id_map();
}

void OffsetBvh::update_face_id_map()
{
    m_bvh_to_mesh_index = bvh_to_mesh_index_map(*m_bvh, m_position_handle, m_pt);
}

double OffsetBvh::sq_dist(const SimpleBVH::VectorMax3d& p, SimpleBVH::VectorMax3d& nearest_point)
    const
{
    double sq_dist;
    m_bvh->nearest_facet(p, nearest_point, sq_dist);
    return sq_dist;
}

double OffsetBvh::sq_dist(const SimpleBVH::VectorMax3d& p) const
{
    SimpleBVH::VectorMax3d nearest_point;
    return sq_dist(p, nearest_point);
}

double OffsetBvh::dist(const SimpleBVH::VectorMax3d& p) const
{
    return std::sqrt(sq_dist(p));
}

std::tuple<double, Tuple> OffsetBvh::sq_dist_and_tuple(
    const SimpleBVH::VectorMax3d& p,
    SimpleBVH::VectorMax3d& nearest_point) const
{
    const Mesh& m = m_position_handle.mesh();

    double sq_dist;
    const int64_t nearest_facet = m_bvh->nearest_facet(p, nearest_point, sq_dist);
    const int64_t fid = m_bvh_to_mesh_index[nearest_facet];

    if (m.top_simplex_type() == m_pt) {
        const Tuple t(0, m.top_simplex_type() == PrimitiveType::Triangle ? 2 : -1, -1, fid);
        return std::make_tuple(sq_dist, t);
    } else {
        const auto faces = m.get_all(m_pt);
        return std::make_tuple(sq_dist, faces[fid]);
    }
}

std::tuple<double, Tuple> OffsetBvh::sq_dist_and_tuple(const SimpleBVH::VectorMax3d& p) const
{
    const Mesh& m = m_position_handle.mesh();

    SimpleBVH::VectorMax3d nearest_point;
    const auto [sq_dist, t] = sq_dist_and_tuple(p, nearest_point);

    return std::make_tuple(sq_dist, t);
}

std::tuple<double, double> OffsetBvh::sq_dist_and_offset_distance(
    const SimpleBVH::VectorMax3d& p,
    SimpleBVH::VectorMax3d& nearest_point) const
{
    const auto [sq_dist, t] = sq_dist_and_tuple(p, nearest_point);

    const Mesh& m = m_f_offset_distance_handle.mesh();
    // const auto acc = m.create_const_accessor<double>(m_f_offset_distance_handle);
    // const double offset_distance = acc.const_scalar_attribute(t);

    // compute offset distance from vertex attribute
    double offset_distance = -1;
    {
        const auto p_acc = m.create_const_accessor<double>(m_position_handle);
        const auto vd_acc = m.create_const_accessor<double>(m_v_offset_distance_handle);
        const auto vertices = m.orient_vertices(t);
        if (vertices.size() == 3) {
            // triangle
            const Eigen::Vector3d p0 = p_acc.const_vector_attribute(vertices[0]);
            const Eigen::Vector3d p1 = p_acc.const_vector_attribute(vertices[1]);
            const Eigen::Vector3d p2 = p_acc.const_vector_attribute(vertices[2]);
            const auto& p = nearest_point;

            const double d0 = vd_acc.const_scalar_attribute(vertices[0]);
            const double d1 = vd_acc.const_scalar_attribute(vertices[1]);
            const double d2 = vd_acc.const_scalar_attribute(vertices[2]);

            const double a0 = 0.5 * (p1 - p).cross(p2 - p).norm();
            const double a1 = 0.5 * (p2 - p).cross(p0 - p).norm();
            const double a2 = 0.5 * (p0 - p).cross(p1 - p).norm();
            const double a = a0 + a1 + a2;
            const double b0 = a0 / a;
            const double b1 = a1 / a;
            const double b2 = a2 / a;

            offset_distance = d0 * b0 + d1 * b1 + d2 * b2;
        } else {
            // edge
            assert(vertices.size() == 2);
            const Eigen::Vector2d p0 = p_acc.const_vector_attribute(vertices[0]);
            const Eigen::Vector2d p1 = p_acc.const_vector_attribute(vertices[1]);
            const auto& p = nearest_point;

            const double d0 = vd_acc.const_scalar_attribute(vertices[0]);
            const double d1 = vd_acc.const_scalar_attribute(vertices[1]);

            const double a0 = (p1 - p).norm();
            const double a1 = (p0 - p).norm();
            const double a = a0 + a1;
            const double b0 = a0 / a;
            const double b1 = a1 / a;

            offset_distance = d0 * b0 + d1 * b1;
        }
    }

    return std::make_tuple(sq_dist, offset_distance);
}

std::tuple<double, double> OffsetBvh::sq_dist_and_offset_distance(
    const SimpleBVH::VectorMax3d& p) const
{
    SimpleBVH::VectorMax3d nearest_point;
    return sq_dist_and_offset_distance(p, nearest_point);
}

std::tuple<double, double, SimpleBVH::VectorMax3d> OffsetBvh::sq_dist_offset_distance_nearest_point(
    const SimpleBVH::VectorMax3d& p) const
{
    SimpleBVH::VectorMax3d nearest_point;
    const auto [sq_dist, offset_distance] = sq_dist_and_offset_distance(p, nearest_point);
    return std::make_tuple(sq_dist, offset_distance, nearest_point);
}

attribute::MeshAttributeHandle OffsetBvh::f_offset_distance_handle() const
{
    return m_f_offset_distance_handle;
}

attribute::MeshAttributeHandle OffsetBvh::v_offset_distance_handle() const
{
    return m_v_offset_distance_handle;
}

void OffsetBvh::update_offset_distance_handle()
{
    Mesh& m = m_position_handle.mesh();
    m_f_offset_distance_handle =
        m.get_attribute_handle<double>(m_f_offset_distance_handle_name, m.top_simplex_type());
    m_v_offset_distance_handle =
        m.get_attribute_handle<double>(m_v_offset_distance_handle_name, PrimitiveType::Vertex);
}

void OffsetBvh::update_position_handle()
{
    Mesh& m = m_position_handle.mesh();
    m_position_handle = m.get_attribute_handle<double>("vertices", PrimitiveType::Vertex);
}

void OffsetBvh::register_offset_distance_attribute(const double offset_distance)
{
    Mesh& m = m_position_handle.mesh();
    m_f_offset_distance_handle = m.register_attribute<double>(
        m_f_offset_distance_handle_name,
        m.top_simplex_type(),
        1,
        false,
        offset_distance);

    m_v_offset_distance_handle = m.register_attribute<double>(
        m_v_offset_distance_handle_name,
        PrimitiveType::Vertex,
        1,
        false,
        offset_distance);
}

const SimpleBVH::BVH& OffsetBvh::bvh() const
{
    return *m_bvh;
}

std::shared_ptr<SimpleBVH::BVH> OffsetBvh::bvh_ptr() const
{
    return m_bvh;
}

Mesh& OffsetBvh::mesh()
{
    return m_position_handle.mesh();
}

const Mesh& OffsetBvh::mesh() const
{
    return m_position_handle.mesh();
}

} // namespace wmtk::components::internal::utils