#include "is_tet_in_offset.hpp"

#include <cmath>
#include <queue>
#include <wmtk/simplex/faces_single_dimension.hpp>
#include <wmtk/utils/bbox_from_mesh.hpp>
#include <wmtk/utils/orient.hpp>

#include "enclosing_sphere.hpp"

namespace wmtk::components::internal::utils {

namespace {
class Cube
{
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor, 3> VecType;

public:
    Cube(const VecType& center, const double length, const int64_t depth = 0)
        : m_center(center)
        , m_length(length)
        , m_depth(depth)
    {}

    Cube(
        const Eigen::MatrixXd& bbox,
        const attribute::MeshAttributeHandle& pos_handle,
        const simplex::Simplex& s)
    {
        const Mesh& m = pos_handle.mesh();

        const auto p_acc = m.create_const_accessor(pos_handle.as<double>());
        const int64_t dim = pos_handle.dimension();


        auto bbox_min = bbox.row(0);
        auto bbox_max = bbox.row(1);

        const double bbox_max_length = (bbox_max - bbox_min).maxCoeff();

        m_center = 0.5 * (bbox_min + bbox_max);
        m_length = (bbox_max - bbox_min).maxCoeff();

        // Eigen::VectorXd tet_bbox_min;
        // Eigen::VectorXd tet_bbox_max;
        //{
        //     const auto vertices = m.orient_vertices(s.tuple());
        //     tet_bbox_min.resize(dim);
        //     tet_bbox_max.resize(dim);
        //     for (int64_t d = 0; d < dim; ++d) {
        //         tet_bbox_min[d] = std::numeric_limits<double>::max();
        //         tet_bbox_max[d] = std::numeric_limits<double>::lowest();
        //     }
        //
        //     for (const Tuple& v : vertices) {
        //         const auto p = p_acc.const_vector_attribute(v);
        //         for (int64_t d = 0; d < dim; ++d) {
        //             tet_bbox_min[d] = std::min(bbox_min[d], p[d]);
        //             tet_bbox_max[d] = std::max(bbox_max[d], p[d]);
        //         }
        //     }
        // }
        // m_center = 0.5 * (tet_bbox_max + tet_bbox_min);
        // m_length = (tet_bbox_max - tet_bbox_min).maxCoeff();
        //
        // Eigen::VectorXd tet_aligned_min;
        //{
        //     tet_aligned_min = bbox_min;
        //     Eigen::VectorXd buffer_max = bbox_max;
        //
        //     for (size_t i = 0; i < 10; ++i) {
        //         Eigen::VectorXd split = 0.5 * (tet_aligned_min + buffer_max);
        //         for (int64_t d = 0; d < dim; ++d) {
        //             if (split[d] < tet_bbox_min[d]) {
        //                 tet_aligned_min[d] = split[d];
        //             } else {
        //                 buffer_max[d] = split[d];
        //             }
        //         }
        //     }
        // }
        // Eigen::VectorXd tet_aligned_max;
        //{
        //     Eigen::VectorXd buffer_min = bbox_min;
        //     tet_aligned_max = bbox_max;
        //
        //     for (size_t i = 0; i < 10; ++i) {
        //         Eigen::VectorXd split = 0.5 * (buffer_min + tet_aligned_max);
        //         for (int64_t d = 0; d < dim; ++d) {
        //             if (split[d] > tet_bbox_max[d]) {
        //                 tet_aligned_max[d] = split[d];
        //             } else {
        //                 buffer_min[d] = split[d];
        //             }
        //         }
        //     }
        // }
        // m_center = 0.5 * (tet_aligned_max + tet_aligned_min);
        // m_length = (tet_aligned_max - tet_aligned_min).maxCoeff();
    }

    /**
     * @brief Radius of the enclosing sphere
     */
    double radius() const
    {
        if (m_center.size() == 2) {
            return 0.5 * m_length * std::sqrt(2);
        } else {
            return 0.5 * m_length * std::sqrt(3);
        }
    }

    int64_t depth() const { return m_depth; }

    double length() const { return m_length; }

    std::vector<Cube> refine() const
    {
        std::vector<Cube> children;
        refine(children);
        return children;
    }

    /**
     * @brief Subdivide cube into 8 new cubes.
     */
    void refine(std::vector<Cube>& cubes) const
    {
        const double l4 = m_length / 4.0;
        const double l2 = m_length / 2.0;
        const int64_t d = m_depth + 1;
        if (m_center.size() == 2) {
            cubes.emplace_back(m_center + l4 * Eigen::Vector2d(1, 1), l2, d);
            cubes.emplace_back(m_center + l4 * Eigen::Vector2d(1, -1), l2, d);
            cubes.emplace_back(m_center + l4 * Eigen::Vector2d(-1, 1), l2, d);
            cubes.emplace_back(m_center + l4 * Eigen::Vector2d(-1, -1), l2, d);
        } else {
            cubes.emplace_back(m_center + l4 * Eigen::Vector3d(1, 1, 1), l2, d);
            cubes.emplace_back(m_center + l4 * Eigen::Vector3d(1, 1, -1), l2, d);
            cubes.emplace_back(m_center + l4 * Eigen::Vector3d(1, -1, 1), l2, d);
            cubes.emplace_back(m_center + l4 * Eigen::Vector3d(1, -1, -1), l2, d);
            cubes.emplace_back(m_center + l4 * Eigen::Vector3d(-1, 1, 1), l2, d);
            cubes.emplace_back(m_center + l4 * Eigen::Vector3d(-1, 1, -1), l2, d);
            cubes.emplace_back(m_center + l4 * Eigen::Vector3d(-1, -1, 1), l2, d);
            cubes.emplace_back(m_center + l4 * Eigen::Vector3d(-1, -1, -1), l2, d);
        }
    }

    void refine(std::queue<Cube>& cubes) const
    {
        const double l4 = m_length / 4.0;
        const double l2 = m_length / 2.0;
        const int64_t d = m_depth + 1;
        if (m_center.size() == 2) {
            cubes.emplace(m_center + l4 * Eigen::Vector2d(1, 1), l2, d);
            cubes.emplace(m_center + l4 * Eigen::Vector2d(1, -1), l2, d);
            cubes.emplace(m_center + l4 * Eigen::Vector2d(-1, 1), l2, d);
            cubes.emplace(m_center + l4 * Eigen::Vector2d(-1, -1), l2, d);
        } else {
            cubes.emplace(m_center + l4 * Eigen::Vector3d(1, 1, 1), l2, d);
            cubes.emplace(m_center + l4 * Eigen::Vector3d(1, 1, -1), l2, d);
            cubes.emplace(m_center + l4 * Eigen::Vector3d(1, -1, 1), l2, d);
            cubes.emplace(m_center + l4 * Eigen::Vector3d(1, -1, -1), l2, d);
            cubes.emplace(m_center + l4 * Eigen::Vector3d(-1, 1, 1), l2, d);
            cubes.emplace(m_center + l4 * Eigen::Vector3d(-1, 1, -1), l2, d);
            cubes.emplace(m_center + l4 * Eigen::Vector3d(-1, -1, 1), l2, d);
            cubes.emplace(m_center + l4 * Eigen::Vector3d(-1, -1, -1), l2, d);
        }
    }

    const VecType& center() const { return m_center; }

private:
    VecType m_center;
    double m_length;
    int64_t m_depth = 0;
};

class TetBVH
{
public:
    TetBVH(const attribute::MeshAttributeHandle& pos_handle, const simplex::Simplex& s)
    {
        constexpr PrimitiveType PV = PrimitiveType::Vertex;
        constexpr PrimitiveType PE = PrimitiveType::Edge;

        const Mesh& m = pos_handle.mesh();

        const auto p_acc = m.create_const_accessor(pos_handle.as<double>());
        const int64_t dim = pos_handle.dimension();

        const PrimitiveType pt_face =
            get_primitive_type_from_id(get_primitive_type_id(s.primitive_type()) - 1);

        const std::vector<Tuple> face_tuples =
            simplex::faces_single_dimension_tuples(m, s, pt_face);

        if (pt_face == PrimitiveType::Triangle) {
            int64_t v_count = 0;
            int64_t f_count = 0;
            assert(p_acc.dimension() == 3);

            V.resize(3 * face_tuples.size(), p_acc.dimension());
            F.resize(face_tuples.size(), 3);

            for (const Tuple& f : face_tuples) {
                auto p0 = p_acc.const_vector_attribute(f);
                auto p1 = p_acc.const_vector_attribute(m.switch_tuple(f, PV));
                auto p2 = p_acc.const_vector_attribute(m.switch_tuples(f, {PE, PV}));

                F.row(f_count) = Eigen::Vector3i(v_count, v_count + 1, v_count + 2);
                V.row(3 * f_count) = p0;
                V.row(3 * f_count + 1) = p1;
                V.row(3 * f_count + 2) = p2;

                v_count += 3;
                ++f_count;
            }
        } else if (pt_face == PrimitiveType::Edge) {
            int64_t v_count = 0;
            int64_t f_count = 0;

            V.resize(2 * face_tuples.size(), p_acc.dimension());
            F.resize(face_tuples.size(), 2);

            for (const Tuple& f : face_tuples) {
                auto p0 = p_acc.const_vector_attribute(f);
                auto p1 = p_acc.const_vector_attribute(m.switch_tuple(f, PV));

                F.row(f_count) = Eigen::Vector2i(v_count, v_count + 1);
                V.row(2 * f_count) = p0;
                V.row(2 * f_count + 1) = p1;

                v_count += 2;
                ++f_count;
            }
        }

        m_bvh = std::make_unique<SimpleBVH::BVH>();
        m_bvh->init(V, F, 1e-10);
    }

    double sq_dist(const SimpleBVH::VectorMax3d& p) const
    {
        double sq_dist;
        SimpleBVH::VectorMax3d nearest_point;
        m_bvh->nearest_facet(p, nearest_point, sq_dist);
        return sq_dist;
    }

private:
    std::unique_ptr<SimpleBVH::BVH> m_bvh;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
};

class TetInside
{
public:
    TetInside(const attribute::MeshAttributeHandle& pos_handle, const simplex::Simplex& s)
    {
        const Mesh& m = pos_handle.mesh();

        const auto p_acc = m.create_const_accessor(pos_handle.as<double>());

        const auto vertices = m.orient_vertices(s.tuple());

        m_pts.reserve(vertices.size());
        for (const Tuple& v : vertices) {
            m_pts.emplace_back(p_acc.const_vector_attribute(v));
        }
    }

    bool is_inside(const SimpleBVH::VectorMax3d& p) const
    {
        using wmtk::utils::wmtk_orient2d;
        using wmtk::utils::wmtk_orient3d;

        if (m_pts.size() == 3) {
            // 2D
            assert(wmtk_orient2d(m_pts[0], m_pts[1], m_pts[2]) == 1);
            bool b0 = wmtk_orient2d(p, m_pts[1], m_pts[2]) == 1;
            bool b1 = wmtk_orient2d(m_pts[0], p, m_pts[2]) == 1;
            bool b2 = wmtk_orient2d(m_pts[0], m_pts[1], p) == 1;
            return b0 && b1 && b2;
        } else {
            assert(m_pts.size() == 4);
            // 3D
            assert(wmtk_orient3d(m_pts[0], m_pts[1], m_pts[2], m_pts[3]) == 1);
            bool b0 = wmtk_orient3d(p, m_pts[1], m_pts[2], m_pts[3]) == 1;
            bool b1 = wmtk_orient3d(m_pts[0], p, m_pts[2], m_pts[3]) == 1;
            bool b2 = wmtk_orient3d(m_pts[0], m_pts[1], p, m_pts[3]) == 1;
            bool b3 = wmtk_orient3d(m_pts[0], m_pts[1], m_pts[2], p) == 1;
            return b0 && b1 && b2 && b3;
        }
    }

private:
    std::vector<SimpleBVH::VectorMax3d> m_pts;
};

} // namespace

bool is_tet_in_offset(
    const attribute::MeshAttributeHandle& pos_handle,
    const SimpleBVH::BVH& bvh,
    const double offset_distance,
    const simplex::Simplex& s)
{
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor, 3> VecType;

    const Mesh& m = pos_handle.mesh();

    const auto p_acc = m.create_const_accessor(pos_handle.as<double>());

    const auto vertices = simplex::faces_single_dimension(m, s, PrimitiveType::Vertex);

    std::vector<VecType> samples;

    VecType mid_point;
    mid_point.setZero(p_acc.dimension(), 1);

    for (const simplex::Simplex& v : vertices) {
        samples.emplace_back(p_acc.const_vector_attribute(v));
        mid_point += samples.back();
    }
    mid_point /= vertices.size();
    samples.emplace_back(mid_point);

    double sq_dist;
    SimpleBVH::VectorMax3d nearest_point;
    for (const auto& sample : samples) {
        bvh.nearest_facet(sample, nearest_point, sq_dist);

        if (sq_dist >= offset_distance * offset_distance) {
            // tet is not fully covered by the offset
            return false;
        }
    }

    return true;
}

bool is_tet_in_offset(
    const attribute::MeshAttributeHandle& pos_handle,
    const OffsetBvh& bvh,
    const simplex::Simplex& s)
{
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor, 3> VecType;

    const Mesh& m = pos_handle.mesh();

    const auto p_acc = m.create_const_accessor(pos_handle.as<double>());

    const auto vertices = m.orient_vertices(s.tuple());
    const auto nv = vertices.size();

    std::vector<VecType> samples;
    samples.reserve(nv + 1);

    VecType mid_point;
    mid_point.setZero(p_acc.dimension(), 1);

    // add vertices and mid point
    for (const Tuple& v : vertices) {
        samples.emplace_back(p_acc.const_vector_attribute(v));
        mid_point += samples.back();
    }
    mid_point /= vertices.size();
    samples.emplace_back(mid_point);

    SimpleBVH::VectorMax3d nearest_point;
    for (const VecType& sample : samples) {
        const auto [sq_dist, offset_distance] = bvh.sq_dist_and_offset_distance(sample);

        if (sq_dist >= offset_distance * offset_distance) {
            // tet is not fully covered by the offset
            return false;
        }
    }

    return true;
}

bool is_tet_in_offset_conservative_sampling(
    const Eigen::MatrixXd& bbox,
    const attribute::MeshAttributeHandle& pos_handle,
    const OffsetBvh& bvh,
    const simplex::Simplex& s,
    const double threshold)
{
    const Mesh& m = pos_handle.mesh();

    const TetBVH tet_bvh(pos_handle, s);
    const TetInside tet_inside(pos_handle, s);

    std::queue<Cube> q;
    q.emplace(bbox, pos_handle, s);

    while (!q.empty()) {
        Cube c = q.front();
        q.pop();

        const double r = c.radius();
        const bool is_inside = tet_inside.is_inside(c.center());

        // check if cube is inside tet or touches it
        if (!is_inside) {
            const double tet_sq_dist = tet_bvh.sq_dist(c.center());
            if (tet_sq_dist > r * r) {
                // cube is not inside tet
                continue;
            }
        }

        const auto [sq_dist, offset_distance] = bvh.sq_dist_and_offset_distance(c.center());
        const double d = std::sqrt(sq_dist);

        // check if cube is outside offset
        if (d - r > offset_distance) {
            return false;
        }

        if (d + r < offset_distance) {
            // cube is inside offset
            // --> continue without refining
            continue;
        }

        // use depth or length for checking if a cube is too small
        if (c.length() < threshold) {
            // cube is too small
            return false;
        }

        // add subdivided cube to queue
        c.refine(q);
    }

    return true;
}

bool is_tet_in_offset_aggressive_sampling(
    const Eigen::MatrixXd& bbox,
    const attribute::MeshAttributeHandle& pos_handle,
    const OffsetBvh& bvh,
    const simplex::Simplex& s,
    const double threshold)
{
    const Mesh& m = pos_handle.mesh();

    const TetBVH tet_bvh(pos_handle, s);
    const TetInside tet_inside(pos_handle, s);

    std::queue<Cube> q;
    q.emplace(bbox, pos_handle, s);

    while (!q.empty()) {
        Cube c = q.front();
        q.pop();

        const double r = c.radius();
        const bool is_inside = tet_inside.is_inside(c.center());

        // check if cube is inside tet or touches it
        if (!is_inside) {
            const double tet_sq_dist = tet_bvh.sq_dist(c.center());
            if (tet_sq_dist > r * r) {
                // cube is not inside tet
                continue;
            }
        }

        const auto [sq_dist, offset_distance] = bvh.sq_dist_and_offset_distance(c.center());
        const double d = std::sqrt(sq_dist);

        // check if cube is fully outside offset
        if (d - r > offset_distance) {
            continue;
        }

        if (d + r < offset_distance) {
            // cube is fully inside offset
            return true;
        }

        // use depth or length for checking if a cube is too small
        if (c.length() < threshold) {
            // cube is too small
            return true;
        }

        // add subdivided cube to queue
        c.refine(q);
    }

    return false;
}

bool is_tet_in_offset_conservative(
    const attribute::MeshAttributeHandle& pos_handle,
    const SimpleBVH::BVH& bvh,
    const double offset_distance,
    const simplex::Simplex& s)
{
    const auto [p, r] = enclosing_sphere(pos_handle, s);

    double sq_dist;
    SimpleBVH::VectorMax3d nearest_point;
    bvh.nearest_facet(p, nearest_point, sq_dist);
    const double dist = std::sqrt(sq_dist);

    return dist + r < offset_distance;
}

} // namespace wmtk::components::internal::utils